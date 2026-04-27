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


# ============================================================
# ==================== USER CONFIG ===========================
# ============================================================

# ---- Runtime / file paths ----
MODEL_CKPT_PATH = "/discover/nobackup/awlee4/GEOSgcm/exp/envs/mlrad/models/mlrad_colfilm1d_qrs.pth"
NORM_PATH       = "/discover/nobackup/awlee4/GEOSgcm/exp/envs/mlrad/norms/NormalizationStats_GEOSgrid_2000-2010.nc"

# ---- Output field name in MAPL export state ----
# EDIT THIS TO MATCH FORTRAN EXPORT SPEC!!!!!!!!!!!!!!!!!!!!!1
EXPORT_NAME = None #"QRS_ML"

# ---- Device ----
# For first integration/debugging, CPU is safer.
DEVICE = "cpu"

# ---- Temporary constant drivers ----
F107_CONST  = 150.0
F107A_CONST = 150.0
KP_CONST    = 3.0
AP_CONST    = 15.0

# ---- Must match training ----
TARGET_VAR = "QRS_TOT"   # or QRS_TOT / QRL_TOT / QJOULE 

HIDDEN_CHANNELS = 256
NUM_BLOCKS      = 5
KERNEL_SIZE     = 5

# ---- Batch size over columns ----
COL_BATCH_SIZE = 256

# ---- Input config; should match training as closely as possible ----
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

# ---- Target transform flags; must match training ----
USE_TARGET_LEV_STANDARDIZE = True
USE_TARGET_ANOMALY_BG      = True

# ---- If LATS/LONS from GEOS are already radians, keep this False ----
GEOS_LATLON_ARE_DEGREES = False


# ============================================================
# ================= GLOBALS / CACHED STATE ===================
# ============================================================

_MODEL = None
_NORM_STATS = None
_Y_MEAN_LEV = None
_Y_STD_LEV = None
_FEATURE_NAMES = None
_COND_DIM = None
_COND_NAMES = None


# ============================================================
# ====================== LOGGING =============================
# ============================================================

def log(msg):
    sys.stderr.write(msg + "\n")
    sys.stderr.flush()


def safe_shape(arr):
    try:
        return str(np.shape(arr))
    except Exception as e:
        return f"<shape failed: {repr(e)}>"


def safe_minmax(arr, name="arr"):
    if arr is None:
        return f"{name}: None"
    try:
        amin = float(np.nanmin(arr))
        amax = float(np.nanmax(arr))
        return f"{name}: shape={np.shape(arr)} min/max={amin:.6g}/{amax:.6g}"
    except Exception as e:
        return f"{name}: shape={safe_shape(arr)} min/max failed: {repr(e)}"


# ============================================================
# =================== NORMALIZATION ==========================
# ============================================================

def load_normalization_stats(nc_path):
    """
    Returns
    -------
    stats : dict
        input normalization stats by feature name
    y_mean_lev : np.ndarray [LM]
    y_std_lev  : np.ndarray [LM]
    """
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
    """
    inputs: (C, LM)
    names : length C
    """
    out = inputs.copy()
    for ch, n in enumerate(names):
        # Harmonics and trig are already bounded features; usually not normalized in your script
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
    """
    pred_trans: (B, LM) or (IM, JM, LM)
    returns same shape in physical units
    """
    out = pred_trans.astype(np.float32)

    if out.ndim == 2:
        # (B, LM)
        if USE_TARGET_LEV_STANDARDIZE:
            out = out * y_std_lev[None, :]
        if USE_TARGET_ANOMALY_BG:
            out = out + y_mean_lev[None, :]

    elif out.ndim == 3:
        # (IM, JM, LM)
        if USE_TARGET_LEV_STANDARDIZE:
            out = out * y_std_lev[None, None, :]
        if USE_TARGET_ANOMALY_BG:
            out = out + y_mean_lev[None, None, :]

    else:
        raise ValueError(f"Unsupported pred_trans.ndim={out.ndim}")

    return out


# ============================================================
# ================= FEATURE / COND HELPERS ===================
# ============================================================

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


# ============================================================
# ====================== MODEL ===============================
# ============================================================

class FiLMResBlock1D(nn.Module):
    def __init__(self, channels, kernel_size=3):
        super().__init__()
        if kernel_size % 2 != 1:
            raise ValueError("kernel_size must be odd")
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
        if kernel_size % 2 != 1:
            raise ValueError("kernel_size must be odd")

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
        # x:    (B, C, LM)
        # cond: (B, D)
        h = self.act(self.in_conv(x))
        gb = self.cond_mlp(cond)  # (B, 2*hidden*num_blocks)
        gb = gb.view(x.size(0), self.num_blocks, 2 * self.hidden_channels)

        for i, block in enumerate(self.blocks):
            gamma, beta = gb[:, i, :].chunk(2, dim=-1)
            gamma = gamma.unsqueeze(-1)  # (B, hidden, 1)
            beta  = beta.unsqueeze(-1)
            h = block(h, gamma, beta)

        out = self.out_conv(h)  # (B, 1, LM)
        return out


def load_model(ckpt_path, in_channels, cond_dim, device=DEVICE):
    if not os.path.exists(ckpt_path):
        raise FileNotFoundError(f"Checkpoint not found: {ckpt_path}")

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

    log(f"[MLRAD] loaded model from {ckpt_path}")
    return model


# ============================================================
# ===================== INIT ONCE ============================
# ============================================================

def init_once():
    global _MODEL, _NORM_STATS, _Y_MEAN_LEV, _Y_STD_LEV, _FEATURE_NAMES, _COND_DIM, _COND_NAMES

    if _MODEL is not None:
        return

    log("[MLRAD] init_once: starting")

    _FEATURE_NAMES = list_selected_feature_names(INPUT_CONFIG)
    _COND_DIM, _COND_NAMES = infer_cond_dim_from_cfg(INPUT_CONFIG)

    log(f"[MLRAD] feature_names ({len(_FEATURE_NAMES)}): {_FEATURE_NAMES}")
    log(f"[MLRAD] cond_names ({_COND_DIM}): {_COND_NAMES}")

    _NORM_STATS, _Y_MEAN_LEV, _Y_STD_LEV = load_normalization_stats(NORM_PATH)
    log(f"[MLRAD] loaded normalization from {NORM_PATH}")
    log(f"[MLRAD] y_mean_lev shape={_Y_MEAN_LEV.shape}, y_std_lev shape={_Y_STD_LEV.shape}")

    _MODEL = load_model(
        MODEL_CKPT_PATH,
        in_channels=len(_FEATURE_NAMES),
        cond_dim=_COND_DIM,
        device=DEVICE,
    )

    log("[MLRAD] init_once: done")


# ============================================================
# ================ GEOS RUNTIME FEATURE BUILDING =============
# ============================================================

def maybe_to_radians(arr2d):
    """
    If GEOS lat/lon are already radians, return as-is.
    If they are degrees, convert to radians.
    """
    arr = np.asarray(arr2d, dtype=np.float32)
    if GEOS_LATLON_ARE_DEGREES:
        return np.deg2rad(arr).astype(np.float32)
    return arr


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
    Returns
    -------
    x_col : np.ndarray, shape (C, LM)
    c_col : np.ndarray, shape (D,)
    """
    LM = T_ijlm.shape[0]

    lat_val = float(lats_rad[i, j])
    lon_val = float(lons_rad[i, j])

    ut_vals  = make_ut_harmonics(hour, cfg["ut_harmonics"])
    doy_vals = make_doy_harmonics(doy,  cfg["doy_harmonics"])

    latlon_sincos = cfg.get("latlon_sincos", False)

    # ------------- x channels -------------
    x_list = []
    for n in feature_names:
        if n in ut_vals:
            x_list.append(np.full((LM,), ut_vals[n], dtype=np.float32))

        elif n in doy_vals:
            x_list.append(np.full((LM,), doy_vals[n], dtype=np.float32))

        elif n == "logp_mid":
            x_list.append(LOGP_ijlm.astype(np.float32))

        elif n == "lat":
            x_list.append(np.full((LM,), lat_val, dtype=np.float32))
        elif n == "lon":
            x_list.append(np.full((LM,), lon_val, dtype=np.float32))

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

    # ------------- conditioning vector -------------
    cond_vals = []
    for n in cond_names:
        if n in ut_vals:
            cond_vals.append(ut_vals[n])
        elif n in doy_vals:
            cond_vals.append(doy_vals[n])

        elif n == "lat":
            cond_vals.append(_normalize_scalar_if_stats(norm_stats, "lat", lat_val))
        elif n == "lon":
            cond_vals.append(_normalize_scalar_if_stats(norm_stats, "lon", lon_val))

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
    f107_const=F107_CONST,
    f107a_const=F107A_CONST,
    kp_const=KP_CONST,
    ap_const=AP_CONST,
    device=DEVICE,
    col_batch_size=COL_BATCH_SIZE,
):
    """
    Assumes arrays come in as:
      T   : (IM, JM, LM)
      PLE : (IM, JM, LM+1)
      LATS/LONS : (IM, JM)

    Returns:
      pred_phys : (IM, JM, LM)
      pred_trans: (IM, JM, LM)
    """
    T = np.asarray(T, dtype=np.float32)
    PLE = np.asarray(PLE, dtype=np.float32)
    LATS = np.asarray(LATS, dtype=np.float32)
    LONS = np.asarray(LONS, dtype=np.float32)

    IM, JM, LM = T.shape

    if PLE.shape != (IM, JM, LM + 1):
        raise RuntimeError(f"Expected PLE shape {(IM, JM, LM+1)}, got {PLE.shape}")
    if LATS.shape != (IM, JM):
        raise RuntimeError(f"Expected LATS shape {(IM, JM)}, got {LATS.shape}")
    if LONS.shape != (IM, JM):
        raise RuntimeError(f"Expected LONS shape {(IM, JM)}, got {LONS.shape}")

    # Convert lat/lon units if needed
    lats_rad = maybe_to_radians(LATS)
    lons_rad = maybe_to_radians(LONS)

    # Mid-level pressure
    P_mid = 0.5 * (PLE[:, :, 0:LM] + PLE[:, :, 1:LM+1])
    P_mid = np.maximum(P_mid, 1.0).astype(np.float32)
    LOGP_MID = np.log(P_mid).astype(np.float32)

    log(f"[MLRAD] {safe_minmax(P_mid, 'P_mid')}")
    log(f"[MLRAD] {safe_minmax(LOGP_MID, 'LOGP_MID')}")

    pred_trans = np.full((IM, JM, LM), np.nan, dtype=np.float32)

    model.eval()
    with torch.no_grad():
        cols_x = []
        cols_c = []
        cols_idx = []

        for i in range(IM):
            for j in range(JM):
                T_col = T[i, j, :]          # (LM,)
                LOGP_col = LOGP_MID[i, j, :]  # (LM,)

                x_col, c_col = build_column_input_and_cond(
                    i=i,
                    j=j,
                    T_ijlm=T_col,
                    LOGP_ijlm=LOGP_col,
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

                if len(cols_x) == col_batch_size or (i == IM - 1 and j == JM - 1):
                    x_batch = np.stack(cols_x, axis=0).astype(np.float32)  # (B, C, LM)
                    c_batch = np.stack(cols_c, axis=0).astype(np.float32)  # (B, D)

                    x_t = torch.from_numpy(x_batch).to(device)
                    c_t = torch.from_numpy(c_batch).to(device)

                    y_t = model(x_t, c_t)                 # (B,1,LM)
                    y_np = y_t.squeeze(1).cpu().numpy()   # (B,LM)

                    for n, (ii, jj) in enumerate(cols_idx):
                        pred_trans[ii, jj, :] = y_np[n, :]

                    cols_x.clear()
                    cols_c.clear()
                    cols_idx.clear()

    pred_phys = invert_target_transform(pred_trans, y_mean_lev, y_std_lev)

    log(f"[MLRAD] {safe_minmax(pred_trans, 'pred_trans')}")
    log(f"[MLRAD] {safe_minmax(pred_phys, 'pred_phys')}")

    return pred_phys, pred_trans


# ============================================================
# ===================== DRIVER CLASS =========================
# ============================================================

class MLRadDriver(UserCode):
    def __init__(self):
        pass


    def init(self, grid_comp, import_state, export_state):
        log("[MLRAD] init (driver)")
        log(f"[MLRAD] sys.version={sys.version}")
        log(f"[MLRAD] PYTHONPATH={os.environ.get('PYTHONPATH')}")
        init_once()


    def run(self, grid_comp, import_state, export_state):
        log(f"[MLRAD] run w/o internal (driver) pid={os.getpid()} t={time.time()}")


    def run_with_internal(self, grid_comp, import_state, export_state, internal_state):
        try:
            log(f"[MLRAD] run_with_internal (driver) pid={os.getpid()} t={time.time()}")

            if internal_state is None:
                raise RuntimeError("[MLRAD] internal_state is None. Use gcrun_with_internal in Fortran.")

            init_once()

            mp = get_MAPLPy()
            IM, JM, LM = mp.grid_dims
            log(f"[MLRAD] IM JM LM = {IM} {JM} {LM}")

            # ---------------- IMPORT ----------------
            T = mp.get_pointer(name="T", state=import_state, dims=[IM, JM, LM])
            PLE = mp.get_pointer(name="PLE", state=import_state, dims=[IM, JM, LM + 1])

            if T is None:
                raise RuntimeError("[MLRAD] T pointer is None")
            if PLE is None:
                raise RuntimeError("[MLRAD] PLE pointer is None")

            # ---------------- INTERNAL ----------------
            LATS = mp.get_pointer(name="MLRAD_LATS", state=internal_state, dims=[IM, JM])
            LONS = mp.get_pointer(name="MLRAD_LONS", state=internal_state, dims=[IM, JM])
            DOY_2D = mp.get_pointer(name="MLRAD_DOY", state=internal_state, dims=[IM, JM])
            HH_2D  = mp.get_pointer(name="MLRAD_HH",  state=internal_state, dims=[IM, JM])

            if LATS is None:
                raise RuntimeError("[MLRAD] MLRAD_LATS pointer is None")
            if LONS is None:
                raise RuntimeError("[MLRAD] MLRAD_LONS pointer is None")
            if DOY_2D is None:
                raise RuntimeError("[MLRAD] MLRAD_DOY pointer is None")
            if HH_2D is None:
                raise RuntimeError("[MLRAD] MLRAD_HH pointer is None")

            DOY = int(DOY_2D[0, 0])
            HH  = int(HH_2D[0, 0])

            log(f"[MLRAD] {safe_minmax(T, 'T')}")
            log(f"[MLRAD] {safe_minmax(PLE, 'PLE')}")
            log(f"[MLRAD] {safe_minmax(LATS, 'LATS')}")
            log(f"[MLRAD] {safe_minmax(LONS, 'LONS')}")
            log(f"[MLRAD] DOY={DOY} HH={HH}")

            # ---------------- INFERENCE ----------------
            pred_phys, pred_trans = run_inference_geos_tile(
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
                f107_const=F107_CONST,
                f107a_const=F107A_CONST,
                kp_const=KP_CONST,
                ap_const=AP_CONST,
                device=DEVICE,
                col_batch_size=COL_BATCH_SIZE,
            )

            # ---------------- WRITE TO EXPORT ----------------
            # NOTE:
            # This assumes the export variable has shape (IM, JM, LM)
            # and expects the same vertical order as pred_phys.
            if EXPORT_NAME is not None:
                out = mp.get_pointer(name=EXPORT_NAME, state=export_state, dims=[IM, JM, LM])
            else:
                out = None

            if out is None:
                log(f"[MLRAD] WARNING: export pointer '{EXPORT_NAME}' is None; skipping write-back.")
            else:
                # If GEOS expects K/s instead of K/day, convert here:
                # out[:, :, :] = pred_phys / 86400.0
                out[:, :, :] = pred_phys 
                log(f"[MLRAD] wrote prediction to export field '{EXPORT_NAME}'")
                log(f"[MLRAD] {safe_minmax(out, EXPORT_NAME)}")

            return

        except Exception as e:
            log(f"[MLRAD] EXCEPTION in run_with_internal: {repr(e)}")
            log(traceback.format_exc())
            raise

    def finalize(self, grid_comp, import_state, export_state):
        log("[MLRAD] finalize (driver)")


CODE = MLRadDriver()
