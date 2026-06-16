# MAPL driver for ML-based radiation heating in GEOS.

import os
import sys
import time
import math
import traceback
from pathlib import Path

from MAPL_PythonBridge import UserCode, get_MAPLPy

import numpy as np
import torch
import torch.nn as nn


# -----------------------------------------------------------------------------
# Runtime config
# -----------------------------------------------------------------------------
TORCH_NUM_THREADS = 1
TORCH_INTEROP_THREADS = 1
COL_BATCH_SIZE = 1024
DEVICE = "cpu"

torch.set_num_threads(TORCH_NUM_THREADS)
try:
    torch.set_num_interop_threads(TORCH_INTEROP_THREADS)
except RuntimeError:
    pass


# -----------------------------------------------------------------------------
# Logging config
# -----------------------------------------------------------------------------
COMPACT_DIAGNOSTICS = True
PRINT_OUTPUT_MINMAX = True


# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
THIS_FILE = Path(__file__).resolve()
MLRAD_COMPONENT_DIR = THIS_FILE.parent.parent
GEOS_SOLAR_DIR = MLRAD_COMPONENT_DIR.parent
GEOS_RADIATION_DIR = GEOS_SOLAR_DIR.parent
GEOS_PHYSICS_DIR = GEOS_RADIATION_DIR.parent
GEOS_AGCM_DIR = GEOS_PHYSICS_DIR.parent

MODEL_DIR = MLRAD_COMPONENT_DIR / "models"
NORM_DIR = MLRAD_COMPONENT_DIR / "norms"

NORM_PATH = NORM_DIR / "NormalizationStats_GEOSgrid_2000-2010.npz"

F107_AP_PATH = (
    GEOS_AGCM_DIR
    / "GEOSsuperdyn_GridComp/@FVdycoreCubed_GridComp/@fvdycore/NRL_MSIS/F107_ap_appended.txt"
)


# -----------------------------------------------------------------------------
# ML field config
# -----------------------------------------------------------------------------
MLRAD_FIELDS = {
    "sw": {
        "enabled": True,
        "export_name": "MLRADSW",
        "ckpt_path": str(MODEL_DIR / "mlrad_colfilm1d_qrs.pth"),
    },
    "lw": {
        "enabled": True,
        "export_name": "MLRADLW",
        "ckpt_path": str(MODEL_DIR / "mlrad_colfilm1d_qrl.pth"),
    },
    "joule": {
        "enabled": True,
        "export_name": "MLRADJH",
        "ckpt_path": str(MODEL_DIR / "mlrad_colfilm1d_qjoule.pth"),
    },
}


# -----------------------------------------------------------------------------
# Space-weather lookup config
# -----------------------------------------------------------------------------
AP_LEVELS = np.array(
    [
        0, 2, 3, 4, 5, 6, 7, 9, 12, 15, 18, 22, 27, 32,
        39, 48, 56, 67, 80, 94, 111, 132, 154, 179,
        207, 236, 300, 400,
    ],
    dtype=np.float32,
)

KP_LEVELS = np.array(
    [
        0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0,
        4.0 / 3.0, 5.0 / 3.0, 2.0,
        7.0 / 3.0, 8.0 / 3.0, 3.0,
        10.0 / 3.0, 11.0 / 3.0, 4.0,
        13.0 / 3.0, 14.0 / 3.0, 5.0,
        16.0 / 3.0, 17.0 / 3.0, 6.0,
        19.0 / 3.0, 20.0 / 3.0, 7.0,
        22.0 / 3.0, 23.0 / 3.0, 8.0,
        25.0 / 3.0, 26.0 / 3.0, 9.0,
    ],
    dtype=np.float32,
)


# -----------------------------------------------------------------------------
# Model architecture config
# -----------------------------------------------------------------------------
HIDDEN_CHANNELS = 256
NUM_BLOCKS = 5
KERNEL_SIZE = 5

P_TRAIN_HPA = np.array(
    [
        9.0000000e-10, 8.5000000e-09, 3.0000000e-07,
        4.9276000e-06, 4.3568900e-05, 2.3067500e-04,
        8.3135800e-04, 2.3185600e-03, 5.2675200e-03,
        1.0000000e-02, 1.6461494e-02, 2.3988098e-02,
        3.3373925e-02, 4.5180100e-02, 5.9846950e-02,
        7.7937295e-02, 1.0010172e-01, 1.2710444e-01,
        1.5973860e-01, 1.9886346e-01, 2.4538486e-01,
        3.0024769e-01, 3.6447639e-01, 4.3915803e-01,
        5.2545737e-01, 6.2453836e-01, 7.3764708e-01,
        8.6607943e-01, 1.0111680e00,
    ],
    dtype=np.float32,
)

P_TRAIN_PA = P_TRAIN_HPA * 100.0
N_ML_LEV = len(P_TRAIN_HPA)

# This matches your current checkpoints:
#   in_channels = 16
#   cond_dim    = 15
#
# Important:
#   training used lev_log = log(lev_geos), where lev_geos is hPa.
#   therefore runtime uses lev_log = log(P_TRAIN_HPA).
INPUT_CONFIG = {
    "ut_harmonics": (1, 2),
    "doy_harmonics": (1, 2),

    "use_lev_log": True,

    "use_lat": True,
    "use_lon": True,
    "latlon_sincos": True,

    "use_f107": True,
    "use_f107a": True,
    "use_kp": True,
    "use_ap": False,

    "use_T": False,
}

USE_TARGET_LEV_STANDARDIZE = True
USE_TARGET_ANOMALY_BG = True


# -----------------------------------------------------------------------------
# Cached runtime objects
# -----------------------------------------------------------------------------
_MODELS = {}
_NORM_STATS = None
_Y_MEAN_LEV = None
_Y_STD_LEV = None
_FEATURE_NAMES = None
_COND_DIM = None
_COND_NAMES = None
_INDEX_TABLE = None

_MLRAD_RANK = int(
    os.environ.get(
        "OMPI_COMM_WORLD_RANK",
        os.environ.get("PMI_RANK", os.environ.get("SLURM_PROCID", "0")),
    )
)


# -----------------------------------------------------------------------------
# Logging helpers
# -----------------------------------------------------------------------------
def log(msg):
    sys.stderr.write(str(msg) + "\n")
    sys.stderr.flush()


def rank0_log(msg):
    if _MLRAD_RANK == 0:
        log(msg)


def arr_minmax_str(arr, name):
    try:
        a = np.asarray(arr)
        return (
            f"{name}: shape={a.shape} "
            f"finite={np.isfinite(a).sum()}/{a.size} "
            f"min={float(np.nanmin(a)):.8e} "
            f"max={float(np.nanmax(a)):.8e} "
            f"mean={float(np.nanmean(a)):.8e}"
        )
    except Exception as e:
        return f"{name}: failed to summarize: {repr(e)}"


# -----------------------------------------------------------------------------
# Space-weather utilities
# -----------------------------------------------------------------------------
def ap_to_kp(ap):
    idx = int(np.argmin(np.abs(AP_LEVELS - np.float32(ap))))
    return float(KP_LEVELS[idx])


def load_f107_ap_table():
    global _INDEX_TABLE

    if _INDEX_TABLE is not None:
        return _INDEX_TABLE

    if not F107_AP_PATH.exists():
        raise FileNotFoundError(f"F107/AP index file not found: {F107_AP_PATH}")

    table = {}
    n_valid = 0

    with F107_AP_PATH.open("r") as f:
        for line in f:
            parts = line.split()

            if not parts or parts[0].startswith("#"):
                continue

            try:
                if len(parts) == 6:
                    year = int(parts[0])
                    doy = int(parts[1])
                    hour = int(parts[2])
                    ap = float(parts[3])
                    f107 = float(parts[4])
                    f107a = float(parts[5])
                elif len(parts) >= 7:
                    year = int(parts[1])
                    doy = int(parts[2])
                    hour = int(parts[3])
                    ap = float(parts[4])
                    f107 = float(parts[5])
                    f107a = float(parts[6])
                else:
                    continue
            except ValueError:
                continue

            key = (year * 1000 + doy) * 24 + hour
            table[key] = (f107, f107a, ap, ap_to_kp(ap))
            n_valid += 1

    if not table:
        raise RuntimeError(f"No valid F107/AP records found in {F107_AP_PATH}")

    _INDEX_TABLE = table
    rank0_log(f"[MLRAD] loaded {n_valid} F107/AP records from {F107_AP_PATH}")

    return _INDEX_TABLE


def get_space_weather_indices(year, doy, hour):
    table = load_f107_ap_table()

    year = int(year)
    doy = int(doy)
    hour = int(round(float(hour)))
    hour = max(0, min(23, hour))

    key = (year * 1000 + doy) * 24 + hour

    if key not in table:
        raise RuntimeError(f"No F107/AP record for year={year}, doy={doy}, hour={hour}")

    f107, f107a, ap, kp = table[key]

    return {
        "f107": f107,
        "f107a": f107a,
        "ap": ap,
        "kp": kp,
    }


# -----------------------------------------------------------------------------
# Normalization
# -----------------------------------------------------------------------------
def load_normalization_stats(npz_path):
    npz_path = Path(npz_path)

    if not npz_path.exists():
        raise FileNotFoundError(f"NORM_PATH not found: {npz_path}")

    data = np.load(str(npz_path), allow_pickle=True)

    stats = {}
    y_mean_lev = None
    y_std_lev = None

    if "var_names" in data and "mean" in data and "std" in data:
        names_raw = data["var_names"]
        means = data["mean"]
        stds = data["std"]

        names = [str(v) for v in names_raw]

        for i, n in enumerate(names):
            stats[n] = {
                "mean": float(means[i]),
                "std": float(stds[i]),
            }

    if "y_mean_lev" in data:
        y_mean_lev = np.array(data["y_mean_lev"], dtype=np.float32)

    if "y_std_lev" in data:
        y_std_lev = np.array(data["y_std_lev"], dtype=np.float32)

    if y_mean_lev is None or y_std_lev is None:
        raise RuntimeError("Normalization file missing y_mean_lev and/or y_std_lev")

    y_std_lev = np.maximum(y_std_lev, 1.0e-6).astype(np.float32)

    return stats, y_mean_lev, y_std_lev


def normalize_scalar_if_stats(stats, name, value):
    if name in stats:
        m = float(stats[name]["mean"])
        s = float(stats[name]["std"]) + 1.0e-6
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
            s = float(stats[n]["std"]) + 1.0e-6
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


# -----------------------------------------------------------------------------
# Feature construction
# -----------------------------------------------------------------------------
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

    if cfg.get("use_lev_log", False):
        names.append("lev_log")

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


def sanity_check_latlon_radians(LATS, LONS):
    """Warn if incoming lat/lon do not look like radians.

    GEOS MLRAD_LATS/MLRAD_LONS are assumed to be radians.
    This function does not convert anything.
    """
    if _MLRAD_RANK != 0:
        return

    lat = np.asarray(LATS, dtype=np.float32)
    lon = np.asarray(LONS, dtype=np.float32)

    lat_abs_max = float(np.nanmax(np.abs(lat)))
    lon_abs_max = float(np.nanmax(np.abs(lon)))

    lat_ok = lat_abs_max <= (math.pi / 2.0 + 1.0e-3)
    lon_ok = lon_abs_max <= (2.0 * math.pi + 1.0e-3)

    if not lat_ok or not lon_ok:
        rank0_log(
            "[MLRAD/WARN] MLRAD_LATS/MLRAD_LONS do not look like radians: "
            f"lat_abs_max={lat_abs_max:.6g}, lon_abs_max={lon_abs_max:.6g}. "
            "This driver assumes radians and will not convert them."
        )


def build_column_input_and_cond(
    i,
    j,
    LEV_LOG_1D,
    T_ML_1D,
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
    LM = LEV_LOG_1D.shape[0]

    if T_ML_1D is not None and T_ML_1D.shape[0] != LM:
        raise RuntimeError(
            f"T_ML_1D length mismatch: len(T_ML_1D)={T_ML_1D.shape[0]}, LM={LM}"
        )

    lat_val = float(lats_rad[i, j])
    lon_val = float(lons_rad[i, j])

    ut_vals = make_ut_harmonics(hour, cfg["ut_harmonics"])
    doy_vals = make_doy_harmonics(doy, cfg["doy_harmonics"])

    x_list = []

    for n in feature_names:
        if n in ut_vals:
            x_list.append(np.full((LM,), ut_vals[n], dtype=np.float32))

        elif n in doy_vals:
            x_list.append(np.full((LM,), doy_vals[n], dtype=np.float32))

        elif n == "lev_log":
            x_list.append(LEV_LOG_1D.astype(np.float32))

        elif n == "T":
            if T_ML_1D is None:
                raise RuntimeError("Feature list contains T, but T_ML_1D is None")
            x_list.append(T_ML_1D.astype(np.float32))

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
            cond_vals.append(normalize_scalar_if_stats(norm_stats, "f107", f107))

        elif n == "f107a":
            cond_vals.append(normalize_scalar_if_stats(norm_stats, "f107a", f107a))

        elif n == "kp":
            cond_vals.append(normalize_scalar_if_stats(norm_stats, "kp", kp))

        elif n == "ap":
            cond_vals.append(normalize_scalar_if_stats(norm_stats, "ap", ap))

        else:
            raise KeyError(f"Unknown feature name in conditioning vector: {n}")

    c_col = np.array(cond_vals, dtype=np.float32)

    return x_col, c_col


# -----------------------------------------------------------------------------
# Model definition
# -----------------------------------------------------------------------------
class FiLMResBlock1D(nn.Module):
    def __init__(self, channels, kernel_size=3):
        super().__init__()

        padding = kernel_size // 2
        self.conv1 = nn.Conv1d(channels, channels, kernel_size=kernel_size, padding=padding)
        self.conv2 = nn.Conv1d(channels, channels, kernel_size=kernel_size, padding=padding)
        self.act = nn.ReLU(inplace=True)

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

        self.in_conv = nn.Conv1d(
            in_channels,
            hidden_channels,
            kernel_size=kernel_size,
            padding=padding,
        )

        self.blocks = nn.ModuleList(
            [
                FiLMResBlock1D(hidden_channels, kernel_size=kernel_size)
                for _ in range(num_blocks)
            ]
        )

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
            beta = beta.unsqueeze(-1)
            h = block(h, gamma, beta)

        return self.out_conv(h)


def load_model(ckpt_path, in_channels, cond_dim, device=DEVICE):
    ckpt_path = str(ckpt_path)

    if not os.path.exists(ckpt_path):
        raise FileNotFoundError(f"Checkpoint not found: {ckpt_path}")

    ckpt = torch.load(ckpt_path, map_location=device)
    sd = ckpt["model"] if "model" in ckpt else ckpt

    if "in_conv.weight" in sd:
        w = sd["in_conv.weight"]
        ckpt_hidden = int(w.shape[0])
        ckpt_in_channels = int(w.shape[1])
        ckpt_kernel = int(w.shape[2])

        rank0_log(
            f"[MLRAD] checkpoint={ckpt_path} "
            f"in_conv.weight={tuple(w.shape)}"
        )

        if ckpt_in_channels != in_channels:
            raise RuntimeError(
                f"Input-channel mismatch for checkpoint {ckpt_path}: "
                f"checkpoint expects {ckpt_in_channels}, "
                f"driver builds {in_channels}. "
                f"Feature names: {_FEATURE_NAMES}"
            )

        if ckpt_hidden != HIDDEN_CHANNELS:
            raise RuntimeError(
                f"HIDDEN_CHANNELS mismatch for checkpoint {ckpt_path}: "
                f"checkpoint has {ckpt_hidden}, driver uses {HIDDEN_CHANNELS}"
            )

        if ckpt_kernel != KERNEL_SIZE:
            raise RuntimeError(
                f"KERNEL_SIZE mismatch for checkpoint {ckpt_path}: "
                f"checkpoint has {ckpt_kernel}, driver uses {KERNEL_SIZE}"
            )

    if "cond_mlp.0.weight" in sd:
        w = sd["cond_mlp.0.weight"]
        ckpt_cond_dim = int(w.shape[1])

        rank0_log(
            f"[MLRAD] checkpoint={ckpt_path} "
            f"cond_mlp.0.weight={tuple(w.shape)}"
        )

        if ckpt_cond_dim != cond_dim:
            raise RuntimeError(
                f"Conditioning-vector mismatch for checkpoint {ckpt_path}: "
                f"checkpoint expects {ckpt_cond_dim}, "
                f"driver builds {cond_dim}. "
                f"Conditioning names: {_COND_NAMES}"
            )

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


# -----------------------------------------------------------------------------
# Initialization
# -----------------------------------------------------------------------------
def init_once():
    global _MODELS, _NORM_STATS, _Y_MEAN_LEV, _Y_STD_LEV
    global _FEATURE_NAMES, _COND_DIM, _COND_NAMES

    if _MODELS:
        return

    rank0_log("[MLRAD] init_once: starting")

    _FEATURE_NAMES = list_selected_feature_names(INPUT_CONFIG)
    _COND_DIM, _COND_NAMES = infer_cond_dim_from_cfg(INPUT_CONFIG)

    _NORM_STATS, _Y_MEAN_LEV, _Y_STD_LEV = load_normalization_stats(NORM_PATH)
    load_f107_ap_table()

    if COMPACT_DIAGNOSTICS:
        rank0_log(f"[MLRAD] NORM_PATH = {NORM_PATH}")
        rank0_log(f"[MLRAD] MODEL_DIR = {MODEL_DIR}")
        rank0_log(f"[MLRAD] feature_names ({len(_FEATURE_NAMES)}) = {_FEATURE_NAMES}")
        rank0_log(f"[MLRAD] cond_names ({_COND_DIM}) = {_COND_NAMES}")

        if "lev_log" in _NORM_STATS:
            rank0_log(
                "[MLRAD] norm[lev_log] "
                f"mean={_NORM_STATS['lev_log']['mean']:.8e} "
                f"std={_NORM_STATS['lev_log']['std']:.8e}"
            )
        else:
            rank0_log("[MLRAD/WARN] no normalization stats for lev_log")

        for name in ("f107", "f107a", "kp"):
            if name in _NORM_STATS:
                rank0_log(
                    f"[MLRAD] norm[{name}] "
                    f"mean={_NORM_STATS[name]['mean']:.8e} "
                    f"std={_NORM_STATS[name]['std']:.8e}"
                )
            else:
                rank0_log(f"[MLRAD/WARN] no normalization stats for {name}")

        rank0_log(arr_minmax_str(_Y_MEAN_LEV, "y_mean_lev [K/day]"))
        rank0_log(arr_minmax_str(_Y_STD_LEV, "y_std_lev [K/day scale]"))
        rank0_log(
            "[MLRAD] lev_log uses log(P_TRAIN_HPA); "
            "log(P_TRAIN_PA) is not used."
        )

        if INPUT_CONFIG.get("use_T", False):
            if "T" in _NORM_STATS:
                rank0_log(
                    "[MLRAD] norm[T] "
                    f"mean={_NORM_STATS['T']['mean']:.8e} "
                    f"std={_NORM_STATS['T']['std']:.8e}"
                )
            else:
                raise RuntimeError(
                    "use_T=True but normalization stats do not contain 'T'. "
                    "Check NormalizationStats_GEOSgrid_2000-2010.npz."
                )

    expected_in_channels = 17 if INPUT_CONFIG.get("use_T", False) else 16
    
    if len(_FEATURE_NAMES) != expected_in_channels:
        raise RuntimeError(
            f"Expected {expected_in_channels} input channels for "
            f"use_T={INPUT_CONFIG.get('use_T', False)}, "
            f"got {len(_FEATURE_NAMES)}: {_FEATURE_NAMES}"
        )

    if _COND_DIM != 15:
        raise RuntimeError(
            f"Expected cond_dim=15 for current config, got {_COND_DIM}: {_COND_NAMES}"
        )

    if len(_Y_MEAN_LEV) != N_ML_LEV or len(_Y_STD_LEV) != N_ML_LEV:
        raise RuntimeError(
            f"Normalization target levels do not match ML levels: "
            f"len(y_mean_lev)={len(_Y_MEAN_LEV)} "
            f"len(y_std_lev)={len(_Y_STD_LEV)} "
            f"N_ML_LEV={N_ML_LEV}"
        )

    for field_key, field_cfg in MLRAD_FIELDS.items():
        if not field_cfg.get("enabled", False):
            continue

        ckpt_path = field_cfg["ckpt_path"]
        export_name = field_cfg["export_name"]

        rank0_log(f"[MLRAD] loading {field_key} model for {export_name}: {ckpt_path}")

        _MODELS[field_key] = load_model(
            ckpt_path,
            in_channels=len(_FEATURE_NAMES),
            cond_dim=_COND_DIM,
            device=DEVICE,
        )

    if not _MODELS:
        raise RuntimeError("No enabled ML radiation fields in MLRAD_FIELDS")

    rank0_log("[MLRAD] init_once: done")


# -----------------------------------------------------------------------------
# Vertical remapping
# -----------------------------------------------------------------------------
def interp_logp_clamp_1d(x_src, y_src, x_tgt):
    """Interpolate in log-pressure space and clamp outside source range.

    x_src:
        source pressure, Pa
    y_src:
        source profile values
    x_tgt:
        target pressure, Pa

    Returns:
        y_src interpolated to x_tgt.
        For target pressures outside source range, use nearest boundary value.
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

    above_top = lxt < lx[0]
    if np.any(above_top):
        out[above_top] = y[0]

    below_bottom = lxt > lx[-1]
    if np.any(below_bottom):
        out[below_bottom] = y[-1]

    return out.astype(np.float32)


def geos_to_ml_temperature(T_col, P_mid_col, p_train_pa):
    """Map GEOS temperature from native GEOS levels to ML pressure levels.

    T_col:
        GEOS temperature column on native GEOS levels.

    P_mid_col:
        GEOS mid-layer pressure column, Pa.

    p_train_pa:
        ML pressure grid, Pa.

    Returns:
        Temperature interpolated to p_train_pa.
    """
    return interp_logp_clamp_1d(
        x_src=P_mid_col,
        y_src=T_col,
        x_tgt=p_train_pa,
    )


def ml_to_geos_heating(heating_ml_col_kday, P_mid_col, p_train_pa, p_ml_bottom_pa):
    """Map ML heating column from ML levels to GEOS levels.

    Input:
        heating_ml_col_kday: K/day on ML pressure grid.

    Output:
        heating_geos_ks: K/s on GEOS levels.
    """
    heating_geos_kday = interp_logp_clamp_1d(
        p_train_pa,
        heating_ml_col_kday,
        P_mid_col,
    )

    heating_geos_kday = np.asarray(heating_geos_kday, dtype=np.float32)

    # Zero below ML lower boundary.
    heating_geos_kday[P_mid_col > p_ml_bottom_pa] = 0.0

    # Training target was converted from K/s to K/day.
    # Convert back to K/s for GEOS.
    heating_geos_ks = heating_geos_kday / 86400.0

    return heating_geos_ks.astype(np.float32)


# -----------------------------------------------------------------------------
# Batch construction and model execution
# -----------------------------------------------------------------------------
def build_ml_tile_inputs(
    T,
    PLE,
    LATS,
    LONS,
    doy,
    hour,
    norm_stats,
    feature_names,
    cond_names,
    f107,
    f107a,
    kp,
    ap,
    device=DEVICE,
    col_batch_size=COL_BATCH_SIZE,
):
    T = np.asarray(T, dtype=np.float32)
    PLE = np.asarray(PLE, dtype=np.float32)
    LATS = np.asarray(LATS, dtype=np.float32)
    LONS = np.asarray(LONS, dtype=np.float32)

    IM, JM, LM = T.shape

    if PLE.shape != (IM, JM, LM + 1):
        raise RuntimeError(
            f"PLE shape mismatch: T.shape={T.shape}, PLE.shape={PLE.shape}"
        )

    P_mid = 0.5 * (PLE[:, :, 0:LM] + PLE[:, :, 1:LM + 1])
    P_mid = np.maximum(P_mid, 1.0e-12).astype(np.float32)

    # GEOS already passes radians. Do not convert.
    lats_rad = LATS
    lons_rad = LONS
    sanity_check_latlon_radians(lats_rad, lons_rad)

    # Match training:
    # lev_geos was in hPa, so lev_log = log(hPa).
    LEV_LOG_1D = np.log(np.maximum(P_TRAIN_HPA, 1.0e-30)).astype(np.float32)

    use_T = INPUT_CONFIG.get("use_T", False)

    batches = []
    cols_x = []
    cols_c = []
    cols_idx = []

    for i in range(IM):
        for j in range(JM):
            if use_T:
                P_mid_col = np.asarray(P_mid[i, j, :], dtype=np.float32)
                T_col = np.asarray(T[i, j, :], dtype=np.float32)

                T_ML_1D = geos_to_ml_temperature(
                    T_col=T_col,
                    P_mid_col=P_mid_col,
                    p_train_pa=P_TRAIN_PA,
                )
            else:
                T_ML_1D = None

            x_col, c_col = build_column_input_and_cond(
                i=i,
                j=j,
                LEV_LOG_1D=LEV_LOG_1D,
                T_ML_1D=T_ML_1D,
                lats_rad=lats_rad,
                lons_rad=lons_rad,
                doy=doy,
                hour=hour,
                f107=f107,
                f107a=f107a,
                kp=kp,
                ap=ap,
                feature_names=feature_names,
                cond_names=cond_names,
                norm_stats=norm_stats,
                cfg=INPUT_CONFIG,
            )

            cols_x.append(x_col)
            cols_c.append(c_col)
            cols_idx.append((i, j))

            if len(cols_x) == col_batch_size or (i == IM - 1 and j == JM - 1):
                x_batch = np.stack(cols_x, axis=0).astype(np.float32)
                c_batch = np.stack(cols_c, axis=0).astype(np.float32)

                x_t = torch.from_numpy(x_batch).to(device)
                c_t = torch.from_numpy(c_batch).to(device)

                batches.append((x_t, c_t, list(cols_idx)))

                cols_x.clear()
                cols_c.clear()
                cols_idx.clear()

    return batches, P_mid


def run_model_from_prebuilt_inputs(
    model,
    batches,
    P_mid,
    y_mean_lev,
    y_std_lev,
    p_ml_bottom_pa,
    device=DEVICE,
):
    IM, JM, LM = P_mid.shape

    pred_trans_ml = np.full((IM, JM, N_ML_LEV), np.nan, dtype=np.float32)

    model.eval()
    inference_context = torch.inference_mode if hasattr(torch, "inference_mode") else torch.no_grad

    with inference_context():
        for x_t, c_t, cols_idx in batches:
            y_t = model(x_t, c_t)
            y_np = y_t.squeeze(1).cpu().numpy()

            for n, (ii, jj) in enumerate(cols_idx):
                pred_trans_ml[ii, jj, :] = y_np[n, :]

    pred_phys_ml_kday = invert_target_transform(pred_trans_ml, y_mean_lev, y_std_lev)

    pred_phys_geos_ks = np.zeros((IM, JM, LM), dtype=np.float32)

    for i in range(IM):
        for j in range(JM):
            pred_phys_geos_ks[i, j, :] = ml_to_geos_heating(
                heating_ml_col_kday=pred_phys_ml_kday[i, j, :],
                P_mid_col=P_mid[i, j, :],
                p_train_pa=P_TRAIN_PA,
                p_ml_bottom_pa=p_ml_bottom_pa,
            )

    return pred_phys_geos_ks


# -----------------------------------------------------------------------------
# MAPL PythonBridge entry point
# -----------------------------------------------------------------------------
class MLRadDriver(UserCode):
    def __init__(self):
        pass

    def init(self, grid_comp, import_state, export_state):
        try:
            rank0_log("[MLRAD] python bridge init")
            init_once()
        except Exception as e:
            log(f"[MLRAD] EXCEPTION in init: {repr(e)}")
            log(traceback.format_exc())
            raise

    def run(self, grid_comp, import_state, export_state):
        pass

    def run_with_internal(self, grid_comp, import_state, export_state, internal_state):
        try:
            t_total0 = time.perf_counter()
            init_once()

            mp = get_MAPLPy()
            IM, JM, LM = mp.grid_dims

            T = mp.get_pointer(name="T", state=import_state, dims=[IM, JM, LM])
            PLE = mp.get_pointer(name="PLE", state=import_state, dims=[IM, JM, LM + 1])

            LATS = mp.get_pointer(name="MLRAD_LATS", state=internal_state, dims=[IM, JM])
            LONS = mp.get_pointer(name="MLRAD_LONS", state=internal_state, dims=[IM, JM])
            YY_2D = mp.get_pointer(name="MLRAD_YY", state=internal_state, dims=[IM, JM])
            DOY_2D = mp.get_pointer(name="MLRAD_DOY", state=internal_state, dims=[IM, JM])
            HH_2D = mp.get_pointer(name="MLRAD_HH", state=internal_state, dims=[IM, JM])
            PBOT_2D = mp.get_pointer(name="MLRAD_PBOT", state=internal_state, dims=[IM, JM])

            YY = int(YY_2D[0, 0])
            DOY = int(DOY_2D[0, 0])
            HH_FLOAT = float(HH_2D[0, 0])
            HH_INDEX = int(round(HH_FLOAT))

            P_ml_bottom_hpa = float(PBOT_2D[0, 0])
            P_ml_bottom_pa = P_ml_bottom_hpa * 100.0

            space_weather = get_space_weather_indices(YY, DOY, HH_INDEX)

            F107_VAL = space_weather["f107"]
            F107A_VAL = space_weather["f107a"]
            KP_VAL = space_weather["kp"]
            AP_VAL = space_weather["ap"]

            t_build0 = time.perf_counter()

            batches, P_mid = build_ml_tile_inputs(
                T=T,
                PLE=PLE,
                LATS=LATS,
                LONS=LONS,
                doy=DOY,
                hour=HH_FLOAT,
                norm_stats=_NORM_STATS,
                feature_names=_FEATURE_NAMES,
                cond_names=_COND_NAMES,
                f107=F107_VAL,
                f107a=F107A_VAL,
                kp=KP_VAL,
                ap=AP_VAL,
                device=DEVICE,
                col_batch_size=COL_BATCH_SIZE,
            )

            t_build = time.perf_counter() - t_build0

            field_timings = []

            for field_key, model in _MODELS.items():
                field_cfg = MLRAD_FIELDS[field_key]
                export_name = field_cfg["export_name"]

                t_field0 = time.perf_counter()

                pred_geos = run_model_from_prebuilt_inputs(
                    model=model,
                    batches=batches,
                    P_mid=P_mid,
                    y_mean_lev=_Y_MEAN_LEV,
                    y_std_lev=_Y_STD_LEV,
                    p_ml_bottom_pa=P_ml_bottom_pa,
                    device=DEVICE,
                )

                t_field = time.perf_counter() - t_field0
                field_timings.append(f"{field_key}={t_field:.3f}s")

                out = mp.get_pointer(name=export_name, state=export_state, dims=[IM, JM, LM])

                if out is None:
                    raise RuntimeError(
                        f"[MLRAD] export pointer is None for {export_name}. "
                        f"Check MAPL_AddExportSpec and runtime allocation."
                    )

                out[:, :, :] = pred_geos[:, :, :]

                if PRINT_OUTPUT_MINMAX and _MLRAD_RANK == 0:
                    rank0_log(arr_minmax_str(out, f"[MLRAD] {export_name} written [K/s]"))

            rank0_log(
                f"[MLRAD_TIMING] rank={_MLRAD_RANK} dims={IM}x{JM}x{LM} "
                f"date={YY}:{DOY}:{HH_FLOAT} "
                f"F107={F107_VAL:.3f} F107A={F107A_VAL:.3f} "
                f"AP={AP_VAL:.3f} KP={KP_VAL:.3f} "
                f"Pbot_hPa={P_ml_bottom_hpa:.6g} "
                f"cols={IM * JM} batches={len(batches)} batch_size={COL_BATCH_SIZE} "
                f"build_inputs={t_build:.3f}s fields={','.join(field_timings)} "
                f"total={time.perf_counter() - t_total0:.3f}s"
            )

        except Exception as e:
            log(f"[MLRAD] EXCEPTION in run_with_internal: {repr(e)}")
            log(traceback.format_exc())
            raise

    def finalize(self, grid_comp, import_state, export_state):
        pass


CODE = MLRadDriver()
