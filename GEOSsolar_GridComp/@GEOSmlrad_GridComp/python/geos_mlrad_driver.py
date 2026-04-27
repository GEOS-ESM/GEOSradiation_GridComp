# MAPL driver for ML-based radiation heating in GEOS.
# Reads GEOS column data, runs a pretrained PyTorch model, and remaps the
# predicted heating tendency back onto GEOS levels for output.

from MAPL_PythonBridge import UserCode, get_MAPLPy
import os
import sys
import time
import math
import traceback

import numpy as np
import netCDF4 as nc
import torch
import torch.nn as nn


MODEL_CKPT_PATH = "/discover/nobackup/awlee4/GEOSgcm/exp/envs/mlrad/models/mlrad_colfilm1d_qrs.pth"
NORM_PATH       = "/discover/nobackup/awlee4/GEOSgcm/exp/envs/mlrad/norms/NormalizationStats_GEOSgrid_2000-2010.nc"

DEVICE = "cpu"

F107_CONST  = 150.0
F107A_CONST = 150.0
KP_CONST    = 3.0
AP_CONST    = 15.0

HIDDEN_CHANNELS = 256
NUM_BLOCKS      = 5
KERNEL_SIZE     = 5
COL_BATCH_SIZE  = 256

P_TRAIN_HPA = np.array([
    9.0000000e-10, 8.5000000e-09, 3.0000000e-07, 4.9276000e-06,
    4.3568900e-05, 2.3067500e-04, 8.3135800e-04, 2.3185600e-03,
    5.2675200e-03, 1.0000000e-02, 1.6461494e-02, 2.3988098e-02,
    3.3373925e-02, 4.5180100e-02, 5.9846950e-02, 7.7937295e-02,
    1.0010172e-01, 1.2710444e-01, 1.5973860e-01, 1.9886346e-01,
    2.4538486e-01, 3.0024769e-01, 3.6447639e-01, 4.3915803e-01,
    5.2545737e-01, 6.2453836e-01, 7.3764708e-01, 8.6607943e-01,
    1.0111680e+00
], dtype=np.float32)

P_TRAIN_PA = P_TRAIN_HPA * 100.0
N_ML_LEV = len(P_TRAIN_PA)

INPUT_CONFIG = {
    "ut_harmonics":  (1, 2),
    "doy_harmonics": (1, 2),
    "use_logp_mid":  True,
    "use_lat":       True,
    "use_lon":       True,
    "latlon_sincos": True,
    "use_f107":      True,
    "use_f107a":     True,
    "use_kp":        True,
    "use_ap":        False,
    "use_T":         True,
}

USE_TARGET_LEV_STANDARDIZE = True
USE_TARGET_ANOMALY_BG      = True
GEOS_LATLON_ARE_DEGREES    = False

# -----------------------------------------------------------------------------
# Module-level cached objects.
# These are loaded once and then reused across driver calls so that the model
# checkpoint and normalization data are not re-read every timestep.
# -----------------------------------------------------------------------------
_MODEL = None
_NORM_STATS = None
_Y_MEAN_LEV = None
_Y_STD_LEV = None
_FEATURE_NAMES = None
_COND_DIM = None
_COND_NAMES = None


def log(msg):
    sys.stderr.write(msg + "\n")
    sys.stderr.flush()


def safe_minmax(arr, name="arr"):
    try:
        return f"{name}: shape={np.shape(arr)} min/max={float(np.nanmin(arr)):.6g}/{float(np.nanmax(arr)):.6g}"
    except Exception as e:
        return f"{name}: failed: {repr(e)}"


def load_normalization_stats(nc_path):
    if not os.path.exists(nc_path):
        raise FileNotFoundError(f"NORM_PATH not found: {nc_path}")

    stats = {}
    y_mean_lev = None
    y_std_lev  = None

    with nc.Dataset(nc_path, "r") as ds:
        if "var_names" in ds.variables and "mean" in ds.variables and "std" in ds.variables:
            names = [str(v) for v in ds.variables["var_names"][:]]
            means = ds.variables["mean"][:]
            stds  = ds.variables["std"][:]
            for i, n in enumerate(names):
                stats[n] = {"mean": float(means[i]), "std": float(stds[i])}

        if "y_mean_lev" in ds.variables:
            y_mean_lev = np.array(ds.variables["y_mean_lev"][:], dtype=np.float32)
        if "y_std_lev" in ds.variables:
            y_std_lev  = np.array(ds.variables["y_std_lev"][:], dtype=np.float32)

    if y_mean_lev is None or y_std_lev is None:
        raise RuntimeError("Normalization file missing y_mean_lev and/or y_std_lev")

    y_std_lev = np.maximum(y_std_lev, 1e-6).astype(np.float32)
    return stats, y_mean_lev, y_std_lev


def _normalize_scalar_if_stats(stats, name, value):
    if name in stats:
        m = float(stats[name]["mean"])
        s = float(stats[name]["std"]) + 1e-6
        return (value - m) / s
    return value


def normalize_inputs_dynamic(inputs, names, stats):
    out = inputs.copy()
    for ch, n in enumerate(names):
        if n.startswith("ut_") or n.startswith("doy_"):
            continue
        if n.endswith("_sin") or n.endswith("_cos"):
            continue
        if n in stats:
            m = float(stats[n]["mean"])
            s = float(stats[n]["std"]) + 1e-6
            out[ch] = (out[ch] - m) / s
    return out


def invert_target_transform(pred_trans, y_mean_lev, y_std_lev):
    out = pred_trans.astype(np.float32)

    if out.ndim == 2:
        if USE_TARGET_LEV_STANDARDIZE:
            out = out * y_std_lev[None, :]
        if USE_TARGET_ANOMALY_BG:
            out = out + y_mean_lev[None, :]
    elif out.ndim == 3:
        if USE_TARGET_LEV_STANDARDIZE:
            out = out * y_std_lev[None, None, :]
        if USE_TARGET_ANOMALY_BG:
            out = out + y_mean_lev[None, None, :]
    else:
        raise ValueError(f"Unsupported pred_trans.ndim={out.ndim}")

    return out


def make_ut_harmonics(hour, ut_harmonics):
    vals = {}
    for h in ut_harmonics:
        ang = (hour / 24.0) * (2.0 * math.pi * h)
        vals[f"ut_cos_h{h}"] = math.cos(ang)
        vals[f"ut_sin_h{h}"] = math.sin(ang)
    return vals


def make_doy_harmonics(doy, doy_harmonics):
    vals = {}
    for h in doy_harmonics:
        ang = (doy / 365.0) * (2.0 * math.pi * h)
        vals[f"doy_cos_h{h}"] = math.cos(ang)
        vals[f"doy_sin_h{h}"] = math.sin(ang)
    return vals


def list_selected_feature_names(cfg=INPUT_CONFIG):
    names = []
    for h in cfg["ut_harmonics"]:
        names += [f"ut_cos_h{h}", f"ut_sin_h{h}"]
    for h in cfg["doy_harmonics"]:
        names += [f"doy_cos_h{h}", f"doy_sin_h{h}"]

    if cfg.get("use_logp_mid", False):
        names.append("logp_mid")

    latlon_sincos = cfg.get("latlon_sincos", False)
    if cfg.get("use_lat", False):
        names += ["lat_sin", "lat_cos"] if latlon_sincos else ["lat"]
    if cfg.get("use_lon", False):
        names += ["lon_sin", "lon_cos"] if latlon_sincos else ["lon"]

    if cfg.get("use_f107", False):
        names.append("f107")
    if cfg.get("use_f107a", False):
        names.append("f107a")
    if cfg.get("use_kp", False):
        names.append("kp")
    if cfg.get("use_ap", False):
        names.append("ap")
    if cfg.get("use_T", False):
        names.append("T")

    return names


def infer_cond_dim_from_cfg(cfg):
    cond_names = []

    for h in cfg.get("ut_harmonics", ()):
        cond_names += [f"ut_cos_h{h}", f"ut_sin_h{h}"]
    for h in cfg.get("doy_harmonics", ()):
        cond_names += [f"doy_cos_h{h}", f"doy_sin_h{h}"]

    latlon_sincos = cfg.get("latlon_sincos", False)
    if cfg.get("use_lat", False):
        cond_names += ["lat_sin", "lat_cos"] if latlon_sincos else ["lat"]
    if cfg.get("use_lon", False):
        cond_names += ["lon_sin", "lon_cos"] if latlon_sincos else ["lon"]

    if cfg.get("use_f107", False):
        cond_names.append("f107")
    if cfg.get("use_f107a", False):
        cond_names.append("f107a")
    if cfg.get("use_kp", False):
        cond_names.append("kp")
    if cfg.get("use_ap", False):
        cond_names.append("ap")

    return len(cond_names), cond_names


class FiLMResBlock1D(nn.Module):
    def __init__(self, channels, kernel_size=3):
        super().__init__()
        padding = kernel_size // 2
        self.conv1 = nn.Conv1d(channels, channels, kernel_size=kernel_size, padding=padding)
        self.conv2 = nn.Conv1d(channels, channels, kernel_size=kernel_size, padding=padding)
        self.act   = nn.ReLU(inplace=True)

    def forward(self, x, gamma, beta):
        out = self.conv1(x)
        out = gamma * out + beta
        out = self.act(out)
        out = self.conv2(out)
        out = out + x
        out = self.act(out)
        return out


class ColumnFiLMCNN1D(nn.Module):
    def __init__(self, in_channels, cond_dim, hidden_channels=64, num_blocks=4, kernel_size=3):
        super().__init__()
        padding = kernel_size // 2
        self.hidden_channels = hidden_channels
        self.num_blocks = num_blocks

        self.in_conv = nn.Conv1d(in_channels, hidden_channels, kernel_size=kernel_size, padding=padding)
        self.blocks = nn.ModuleList([
            FiLMResBlock1D(hidden_channels, kernel_size=kernel_size) for _ in range(num_blocks)
        ])
        self.out_conv = nn.Conv1d(hidden_channels, 1, kernel_size=1)

        self.cond_mlp = nn.Sequential(
            nn.Linear(cond_dim, hidden_channels * 2),
            nn.ReLU(inplace=True),
            nn.Linear(hidden_channels * 2, hidden_channels * 2 * num_blocks),
        )
        self.act = nn.ReLU(inplace=True)

    def forward(self, x, cond):
        h = self.act(self.in_conv(x))
        gb = self.cond_mlp(cond)
        gb = gb.view(x.size(0), self.num_blocks, 2 * self.hidden_channels)

        for i, block in enumerate(self.blocks):
            gamma, beta = gb[:, i, :].chunk(2, dim=-1)
            gamma = gamma.unsqueeze(-1)
            beta  = beta.unsqueeze(-1)
            h = block(h, gamma, beta)

        return self.out_conv(h)


def load_model(ckpt_path, in_channels, cond_dim, device=DEVICE):
    ckpt = torch.load(ckpt_path, map_location=device)
    sd = ckpt["model"] if "model" in ckpt else ckpt

    model = ColumnFiLMCNN1D(
        in_channels=in_channels,
        cond_dim=cond_dim,
        hidden_channels=HIDDEN_CHANNELS,
        num_blocks=NUM_BLOCKS,
        kernel_size=KERNEL_SIZE,
    ).to(device)

    model.load_state_dict(sd, strict=True)
    model.eval()
    return model


def init_once():
    """
    initialize the ML model and normalization metadata.

    MAPL may call the driver multiple times, so we cache these objects at the
    module level and only load them once per process.
    """
    global _MODEL, _NORM_STATS, _Y_MEAN_LEV, _Y_STD_LEV, _FEATURE_NAMES, _COND_DIM, _COND_NAMES

    if _MODEL is not None:
        return

    log("[MLRAD] init_once: starting")
    _FEATURE_NAMES = list_selected_feature_names(INPUT_CONFIG)
    _COND_DIM, _COND_NAMES = infer_cond_dim_from_cfg(INPUT_CONFIG)
    _NORM_STATS, _Y_MEAN_LEV, _Y_STD_LEV = load_normalization_stats(NORM_PATH)
    _MODEL = load_model(
        MODEL_CKPT_PATH,
        in_channels=len(_FEATURE_NAMES),
        cond_dim=_COND_DIM,
        device=DEVICE,
    )
    if len(_Y_MEAN_LEV) != N_ML_LEV or len(_Y_STD_LEV) != N_ML_LEV:
        raise RuntimeError(
            f"Normalization target levels do not match ML levels: "
            f"len(y_mean_lev)={len(_Y_MEAN_LEV)} "
            f"len(y_std_lev)={len(_Y_STD_LEV)} "
            f"N_ML_LEV={N_ML_LEV}"
        )

def maybe_to_radians(arr2d):
    """
    Convert latitude/longitude to radians only if the incoming GEOS fields are
    stored in degrees. The ML feature builder assumes radians for trig features.
    """
    arr = np.asarray(arr2d, dtype=np.float32)
    if GEOS_LATLON_ARE_DEGREES:
        return np.deg2rad(arr).astype(np.float32)
    return arr

def interp_extrap_logp_1d(x_src, y_src, x_tgt):
    """
    Interpolate y(x) in log-pressure space with linear extrapolation on both ends.
    x_src, x_tgt are pressures in Pa. Works best when x_src is monotonic.
    """
    x_src = np.asarray(x_src, dtype=np.float64)
    y_src = np.asarray(y_src, dtype=np.float64)
    x_tgt = np.asarray(x_tgt, dtype=np.float64)

    good = np.isfinite(x_src) & np.isfinite(y_src) & (x_src > 0.0)
    x = x_src[good]
    y = y_src[good]

    if x.size < 2:
        if x.size == 1:
            return np.full_like(x_tgt, y[0], dtype=np.float32)
        return np.full_like(x_tgt, np.nan, dtype=np.float32)

    lx = np.log(x)
    lxt = np.log(np.maximum(x_tgt, 1.0e-30))

    order = np.argsort(lx)
    lx = lx[order]
    y = y[order]

    out = np.interp(lxt, lx, y).astype(np.float64)

    left = lxt < lx[0]
    if np.any(left):
        m = (y[1] - y[0]) / (lx[1] - lx[0])
        out[left] = y[0] + m * (lxt[left] - lx[0])

    right = lxt > lx[-1]
    if np.any(right):
        m = (y[-1] - y[-2]) / (lx[-1] - lx[-2])
        out[right] = y[-1] + m * (lxt[right] - lx[-1])

    return out.astype(np.float32)


def interp_logp_top_clamp_1d(x_src, y_src, x_tgt):
    """
    Interpolate y(x) in log-pressure space.
    Above the GEOS top (lower pressure than min x_src), clamp to top value.
    Below the GEOS bottom (higher pressure than max x_src), clamp to bottom value.
    x_src, x_tgt are pressures in Pa.
    """
    x_src = np.asarray(x_src, dtype=np.float64)
    y_src = np.asarray(y_src, dtype=np.float64)
    x_tgt = np.asarray(x_tgt, dtype=np.float64)

    good = np.isfinite(x_src) & np.isfinite(y_src) & (x_src > 0.0)
    x = x_src[good]
    y = y_src[good]

    if x.size < 2:
        if x.size == 1:
            return np.full_like(x_tgt, y[0], dtype=np.float32)
        return np.full_like(x_tgt, np.nan, dtype=np.float32)

    lx = np.log(x)
    lxt = np.log(np.maximum(x_tgt, 1.0e-30))

    order = np.argsort(lx)
    lx = lx[order]
    y = y[order]

    out = np.interp(lxt, lx, y).astype(np.float64)

    # target pressure smaller than GEOS top pressure => above model top
    above_top = lxt < lx[0]
    if np.any(above_top):
        out[above_top] = y[0]

    # target pressure larger than GEOS bottom pressure => below model bottom
    below_bottom = lxt > lx[-1]
    if np.any(below_bottom):
        out[below_bottom] = y[-1]

    return out.astype(np.float32)


def geos_to_ml_temperature(T_col, P_mid_col, p_train_pa):
    #return interp_extrap_logp_1d(P_mid_col, T_col, p_train_pa)

    # Once you run the true 190-level configuration, 
    # you should no longer need this clamp for those top levels.
    return interp_logp_top_clamp_1d(P_mid_col, T_col, p_train_pa)


def ml_to_geos_heating(heating_ml_col, P_mid_col, p_train_pa, p_ml_bottom_pa):
    """
    Remap ML-predicted heating from the ML pressure grid back to GEOS levels.

    The ML output is defined only over the supported pressure range, so heating
    is set to zero below p_ml_bottom_pa. Units are converted from K/day to K/s
    before returning to GEOS.
    """
    heating_geos = interp_extrap_logp_1d(p_train_pa, heating_ml_col, P_mid_col)
    heating_geos = np.asarray(heating_geos, dtype=np.float32)

    # mask below the allowed ML region
    heating_geos[P_mid_col > p_ml_bottom_pa] = 0.0

    # convert K/day -> K/s
    heating_geos = heating_geos / 86400.0

    return heating_geos


def build_column_input_and_cond(
    i, j,
    T_ijlm,
    LOGP_ijlm,
    lats_rad,
    lons_rad,
    doy,
    hour,
    f107,
    f107a,
    kp,
    ap,
    feature_names,
    cond_names,
    norm_stats,
    cfg,
):
    """
    Build one column's ML input tensor and FiLM conditioning vector.

    x_col contains level-dependent features on the ML pressure grid.
    c_col contains column-wise conditioning features shared across levels
    (time harmonics, location, and space-weather drivers).
    """
    LM = T_ijlm.shape[0]

    lat_val = float(lats_rad[i, j])
    lon_val = float(lons_rad[i, j])

    ut_vals  = make_ut_harmonics(hour, cfg["ut_harmonics"])
    doy_vals = make_doy_harmonics(doy,  cfg["doy_harmonics"])

    x_list = []
    for n in feature_names:
        if n in ut_vals:
            x_list.append(np.full((LM,), ut_vals[n], dtype=np.float32))
        elif n in doy_vals:
            x_list.append(np.full((LM,), doy_vals[n], dtype=np.float32))
        elif n == "logp_mid":
            x_list.append(LOGP_ijlm.astype(np.float32))
        elif n == "lat_sin":
            x_list.append(np.full((LM,), math.sin(lat_val), dtype=np.float32))
        elif n == "lat_cos":
            x_list.append(np.full((LM,), math.cos(lat_val), dtype=np.float32))
        elif n == "lon_sin":
            x_list.append(np.full((LM,), math.sin(lon_val), dtype=np.float32))
        elif n == "lon_cos":
            x_list.append(np.full((LM,), math.cos(lon_val), dtype=np.float32))
        elif n == "f107":
            x_list.append(np.full((LM,), f107, dtype=np.float32))
        elif n == "f107a":
            x_list.append(np.full((LM,), f107a, dtype=np.float32))
        elif n == "kp":
            x_list.append(np.full((LM,), kp, dtype=np.float32))
        elif n == "ap":
            x_list.append(np.full((LM,), ap, dtype=np.float32))
        elif n == "T":
            x_list.append(T_ijlm.astype(np.float32))
        else:
            raise KeyError(f"Unknown feature name in x-channels: {n}")

    x_col = np.stack(x_list, axis=0).astype(np.float32)
    x_col = normalize_inputs_dynamic(x_col, feature_names, norm_stats)

    cond_vals = []
    for n in cond_names:
        if n in ut_vals:
            cond_vals.append(ut_vals[n])
        elif n in doy_vals:
            cond_vals.append(doy_vals[n])
        elif n == "lat_sin":
            cond_vals.append(math.sin(lat_val))
        elif n == "lat_cos":
            cond_vals.append(math.cos(lat_val))
        elif n == "lon_sin":
            cond_vals.append(math.sin(lon_val))
        elif n == "lon_cos":
            cond_vals.append(math.cos(lon_val))
        elif n == "f107":
            cond_vals.append(_normalize_scalar_if_stats(norm_stats, "f107", f107))
        elif n == "f107a":
            cond_vals.append(_normalize_scalar_if_stats(norm_stats, "f107a", f107a))
        elif n == "kp":
            cond_vals.append(_normalize_scalar_if_stats(norm_stats, "kp", kp))
        elif n == "ap":
            cond_vals.append(_normalize_scalar_if_stats(norm_stats, "ap", ap))
        else:
            raise KeyError(f"Unknown feature name in conditioning vector: {n}")

    c_col = np.array(cond_vals, dtype=np.float32)
    return x_col, c_col


def run_inference_geos_tile(
    model,
    T,
    PLE,
    LATS,
    LONS,
    doy,
    hour,
    norm_stats,
    feature_names,
    cond_names,
    y_mean_lev,
    y_std_lev,
    p_ml_bottom_pa,
    f107_const=F107_CONST,
    f107a_const=F107A_CONST,
    kp_const=KP_CONST,
    ap_const=AP_CONST,
    device=DEVICE,
    col_batch_size=COL_BATCH_SIZE,
):
    """
    Run ML inference over one GEOS tile.

    Workflow:
      1. Compute GEOS layer-midpoint pressure.
      2. Remap each column onto the ML pressure grid.
      3. Build feature tensors and conditioning vectors.
      4. Run batched PyTorch inference.
      5. Invert target normalization.
      6. Remap predicted heating back to GEOS levels.
    """
    T = np.asarray(T, dtype=np.float32)
    PLE = np.asarray(PLE, dtype=np.float32)
    LATS = np.asarray(LATS, dtype=np.float32)
    LONS = np.asarray(LONS, dtype=np.float32)

    IM, JM, LM = T.shape
    
    # Convert edge pressure to layer-midpoint pressure for GEOS vertical levels.
    P_mid = 0.5 * (PLE[:, :, 0:LM] + PLE[:, :, 1:LM+1])
    P_mid = np.maximum(P_mid, 1.0e-12).astype(np.float32)

    log(f"[MLRAD] {safe_minmax(P_mid, 'P_mid')}")

    lats_rad = maybe_to_radians(LATS)
    lons_rad = maybe_to_radians(LONS)

    pred_trans_ml = np.full((IM, JM, N_ML_LEV), np.nan, dtype=np.float32)
    pred_phys_ml  = np.full((IM, JM, N_ML_LEV), np.nan, dtype=np.float32)
    pred_phys_geos = np.zeros((IM, JM, LM), dtype=np.float32)
    # pred_trans_ml : raw network output in normalized target space
    # pred_phys_ml  : ML-grid heating after inverse normalization
    # pred_phys_geos: heating remapped back to native GEOS levels

    model.eval()
    with torch.no_grad():
        cols_x = []
        cols_c = []
        cols_idx = []

        for i in range(IM):
            for j in range(JM):
                P_mid_col = np.asarray(P_mid[i, j, :], dtype=np.float32)
                T_col     = np.asarray(T[i, j, :], dtype=np.float32)

                T_ml_col = geos_to_ml_temperature(T_col, P_mid_col, P_TRAIN_PA)
                if i == 0 and j == 0:
                    log(f"[MLRAD] top GEOS P_mid_col[:5] = {P_mid_col[:5]}")
                    log(f"[MLRAD] top train P_TRAIN_PA[:5] = {P_TRAIN_PA[:5]}")
                    log(f"[MLRAD] top GEOS T_col[:5] = {T_col[:5]}")
                    log(f"[MLRAD] top ML   T_ml_col[:5] = {T_ml_col[:5]}")
                LOGP_ml_col = np.log(np.maximum(P_TRAIN_PA, 1.0e-30)).astype(np.float32)

                x_col, c_col = build_column_input_and_cond(
                    i=i,
                    j=j,
                    T_ijlm=T_ml_col,
                    LOGP_ijlm=LOGP_ml_col,
                    lats_rad=lats_rad,
                    lons_rad=lons_rad,
                    doy=doy,
                    hour=hour,
                    f107=f107_const,
                    f107a=f107a_const,
                    kp=kp_const,
                    ap=ap_const,
                    feature_names=feature_names,
                    cond_names=cond_names,
                    norm_stats=norm_stats,
                    cfg=INPUT_CONFIG,
                )

                cols_x.append(x_col)
                cols_c.append(c_col)
                cols_idx.append((i, j))

                # Accumulate columns and run inference in batches to reduce Python overhead
                # and avoid calling the network one column at a time.
                if len(cols_x) == col_batch_size or (i == IM - 1 and j == JM - 1):
                    x_batch = np.stack(cols_x, axis=0).astype(np.float32)   # (B, C, 29)
                    c_batch = np.stack(cols_c, axis=0).astype(np.float32)   # (B, cond)

                    x_t = torch.from_numpy(x_batch).to(device)
                    c_t = torch.from_numpy(c_batch).to(device)

                    y_t = model(x_t, c_t)                                   # (B, 1, 29)
                    y_np = y_t.squeeze(1).cpu().numpy()                     # (B, 29)

                    for n, (ii, jj) in enumerate(cols_idx):
                        pred_trans_ml[ii, jj, :] = y_np[n, :]

                    cols_x.clear()
                    cols_c.clear()
                    cols_idx.clear()

    pred_phys_ml = invert_target_transform(pred_trans_ml, y_mean_lev, y_std_lev)

    for i in range(IM):
        for j in range(JM):
            pred_phys_geos[i, j, :] = ml_to_geos_heating(
                heating_ml_col=pred_phys_ml[i, j, :],
                P_mid_col=P_mid[i, j, :],
                p_train_pa=P_TRAIN_PA,
                p_ml_bottom_pa=p_ml_bottom_pa,
            )

    return pred_phys_geos, pred_phys_ml, pred_trans_ml


class MLRadDriver(UserCode):
    def __init__(self):
        pass

    def init(self, grid_comp, import_state, export_state):
        try:
            log("[MLRAD] Stage3 init (driver)")
            init_once()
        except Exception as e:
            log(f"[MLRAD] EXCEPTION in init: {repr(e)}")
            log(traceback.format_exc())
            raise

    def run(self, grid_comp, import_state, export_state):
        log(f"[MLRAD] Stage3 run w/o internal pid={os.getpid()} t={time.time()}")

    def run_with_internal(self, grid_comp, import_state, export_state, internal_state):
        """
        Main MAPL entry point for this driver.
    
        Reads GEOS state variables and auxiliary internal fields, runs the ML
        heating prediction on the current tile, and writes the result to the
        export-state field MLRADSW.
        """
        try:
            log(f"[MLRAD] Stage3 run_with_internal pid={os.getpid()} t={time.time()}")
            init_once()

            mp = get_MAPLPy()
            IM, JM, LM = mp.grid_dims
            log(f"[MLRAD] IM JM LM = {IM} {JM} {LM}")

            # Import state: prognostic GEOS fields needed by the ML model.
            T = mp.get_pointer(name="T", state=import_state, dims=[IM, JM, LM])
            PLE = mp.get_pointer(name="PLE", state=import_state, dims=[IM, JM, LM + 1])

            # Internal state: precomputed metadata and driver-control fields.
            LATS = mp.get_pointer(name="MLRAD_LATS", state=internal_state, dims=[IM, JM])
            LONS = mp.get_pointer(name="MLRAD_LONS", state=internal_state, dims=[IM, JM])
            DOY_2D = mp.get_pointer(name="MLRAD_DOY", state=internal_state, dims=[IM, JM])
            HH_2D  = mp.get_pointer(name="MLRAD_HH",  state=internal_state, dims=[IM, JM])
            PBOT_2D = mp.get_pointer(name="MLRAD_PBOT", state=internal_state, dims=[IM, JM])

            DOY = int(DOY_2D[0, 0])
            HH  = int(HH_2D[0, 0])
            P_ml_bottom_hpa = float(PBOT_2D[0, 0])
            P_ml_bottom_pa  = P_ml_bottom_hpa * 100.0
            log(f"[MLRAD] P_ml_bottom_hpa={P_ml_bottom_hpa:.6g}  P_ml_bottom_pa={P_ml_bottom_pa:.6g}")

            log(f"[MLRAD] {safe_minmax(T, 'T')}")
            log(f"[MLRAD] {safe_minmax(PLE, 'PLE')}")
            log(f"[MLRAD] {safe_minmax(LATS, 'LATS')}")
            log(f"[MLRAD] {safe_minmax(LONS, 'LONS')}")
            log(f"[MLRAD] DOY={DOY} HH={HH}")

            pred_geos, pred_phys_ml, pred_trans_ml = run_inference_geos_tile(
                model=_MODEL,
                T=T,
                PLE=PLE,
                LATS=LATS,
                LONS=LONS,
                doy=DOY,
                hour=HH,
                norm_stats=_NORM_STATS,
                feature_names=_FEATURE_NAMES,
                cond_names=_COND_NAMES,
                y_mean_lev=_Y_MEAN_LEV,
                y_std_lev=_Y_STD_LEV,
                p_ml_bottom_pa=P_ml_bottom_pa,
                f107_const=F107_CONST,
                f107a_const=F107A_CONST,
                kp_const=KP_CONST,
                ap_const=AP_CONST,
                device=DEVICE,
                col_batch_size=COL_BATCH_SIZE,
            )

            log(f"[MLRAD] {safe_minmax(pred_trans_ml, 'pred_trans_ml')}")
            log(f"[MLRAD] {safe_minmax(pred_phys_ml, 'pred_phys_ml')}")
            log(f"[MLRAD] {safe_minmax(pred_geos, 'pred_geos')}")

            # Write heating tendency back to the GEOS export state.
            MLRADSW = mp.get_pointer(name="MLRADSW", state=export_state, dims=[IM, JM, LM])
            MLRADSW[:, :, :] = pred_geos[:, :, :]

        except Exception as e:
            log(f"[MLRAD] EXCEPTION in run_with_internal: {repr(e)}")
            log(traceback.format_exc())
            raise

    def finalize(self, grid_comp, import_state, export_state):
        log("[MLRAD] finalize (driver)")


CODE = MLRadDriver()
