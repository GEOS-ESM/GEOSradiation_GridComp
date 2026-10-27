# GEOSradiation regression tests

Standalone MAPL3 regression cases for the radiation components (IRRAD,
SOLAR, SATSIM and the RADIATION container). Each case runs one component
(or RADIATION with its children) under a configurable `Root`, with every
import supplied from a restarted "DataFor*" gridcomp, and compares the
resulting checkpoints against a stored baseline.

| Case        | Component(s) run           | States compared by ctest                         |
|-------------|----------------------------|--------------------------------------------------|
| `irrad-sa`  | IRRAD                      | IRRAD import, internal, export                   |
| `solar-sa`  | SOLAR                      | SOLAR import, internal, export                   |
| `satsim-sa` | RADIATION stand-in -> SATSIM (ISCCP only) | SATSIM import, export                            |
| `radiation` | RADIATION -> IRRAD + SOLAR + SATSIM (ISCCP only) | RADIATION import, export; IRRAD and SOLAR import, internal, export; SATSIM import, export |

The cases are listed in `cases.txt` and registered by `CMakeLists.txt`
through `esma_add_regression_tests` (ctest label `REGRESSION`, test name =
case name).

## Common setup

All four cases share:

- **Geometry**: C12 cubed sphere, 181 levels (C12L181), model vertical grid
  defined by the `PLE` field (`root.yaml`).
- **Layout**: 6 MPI tasks (`mapl.yaml`, `run_case.cmake`).
- **Clock**: one 20-minute step (`PT20M`) from `2000-04-14T22:00:00`
  (`cap_restart.yaml`); the final checkpoint is written at 22:20
  (`checkpoints/last`). The `start`/`stop` in `cap_driver.yaml` only bound
  the clock.
- **No ExtData, no History output**: `run_extdata: false`; History is off
  except in `satsim-sa` and `radiation` (see below), where it never writes
  a file.
- **Inputs**: each case's `DataFor*` gridcomp
  (`data-for-*.yaml`, `libconfigurable_gridcomp`) holds every import as an
  internal field, restarted from
  `checkpoints/2000-04-14T22:00:00/DataFor*_internal.nc`, and is connected
  field-by-field to the component under test in `root.yaml`.
- **Radiation schemes** (IRRAD/SOLAR): RRTMGP for both LW and SW
  (`USE_RRTMGP_IRRAD: 1`, `USE_RRTMGP_SORAD: 1`), RRTMGP v1.8 data files
  from `ExtData/g5gcm/radiation/RRTMGP/v1.8/`, `NUM_BANDS: 30`.
- **Aerosols off**: `FakeGOCART_GridComp` (in this directory) provides an
  empty `AERO` state with `implements_aerosol_optics_method = .false.`,
  so IRRAD and SOLAR skip aerosol optics. This matches the MAPL2 baseline
  run's `AERO_PROVIDER: none` (see below).
- **Exports**: `activate_all_exports: true` for IRRAD, SOLAR and
  RADIATION; SATSIM's exports are activated through History instead.

## Cases

### `irrad-sa`

IRRAD alone. `DataForIrrad` supplies 24 imports (`PLE`, `T`, `QV`,
`QL/QI/QR/QS/QG`, `RL/RI/RR/RS/RG`, `O3`, `CH4`, `N2O`, `CFC11`, `CFC12`,
`HCFC22`, `FCLD`, `TS`, `EMIS`, `PREF`, `TSINST`); `AERO` comes from
FakeGOCART. IRRAD's own internal state is restarted from
`checkpoints/2000-04-14T22:00:00/IRRAD_internal.nc`.

Assumptions (`irrad.yaml`):

- `IRRAD_DT: 1200`, equal to the clock step, so the LW alarm rings and
  `LW_Driver` runs on the single step.
- `CO2: -1.0`; the 3-D `CO2` import is not provided (commented out in
  `data-for-irrad.yaml` and `root.yaml`, needed only for `USE_CO2_3D`).
- `RATS_DIAGNOSTICS: []`: no RATS diagnostics.

### `solar-sa`

SOLAR alone. `DataForSolar` supplies 23 imports (as `irrad-sa` minus
`CFC11/CFC12/HCFC22`, `EMIS`, `TSINST`, plus the four surface albedos
`ALBVR/ALBVF/ALBNR/ALBNF`); `AERO` comes from FakeGOCART. SOLAR's internal
state is restarted from `checkpoints/2000-04-14T22:00:00/SOLAR_internal.nc`.

Assumptions (`solar.yaml`, which mirrors the `AGCM.rc` of the MAPL2
baseline run):

- `SOLAR_DT: 1200`, equal to the clock step, so the SW refresh runs on the
  single step.
- Orbit: `EOT: true`, `ORBIT_ANAL2B: true`; solar variability from NRLSSI2
  (`ISOLVAR: 2`, `USE_NRLSSI2: true`, `SOLAR_CYCLE_FILE_NAME:
  ExtData/g5gcm/solar/NRLSSI2.txt`); `SOLAR_CONSTANT: -1`, `CO2: -1`.
- `RAD_CONDENSATE_INHOMOGENEITY: 1` is set here directly; in the full model
  the RADIATION parent sets it.
- `activate_all_exports: true` turns band output on for the supported bands
  08-11, so `[IO]SRB{08..11}RGN` are computed. The MAPL2 baseline run must
  request `[IO]SRB{08..11}RG` from SOLAR in `HISTORY.rc` for these to match.
- With `CALLED_LAST` at its default, the export update runs before the
  refresh and `[IO]SRBbbRG` are computed from the restarted internal state.
  `SOLAR_internal.nc` must therefore come from the same MAPL2 capture as
  the baseline.

### `satsim-sa`

SATSIM, ISCCP simulator only, as the only child of a stand-in RADIATION
(`radiation.yaml`, a `libconfigurable_gridcomp` that only declares the
`USE_SATSIM*` flags, which SATSIM inherits as in the full model). SATSIM's
imports propagate up to RADIATION, where `root.yaml` connects them.
`DataForSatsim` supplies SATSIM's 19 imports (`T`, `QV`, `FCLD`, `RL/RI/RG/RS/RR`, `PLE`, `ZLE`, `QL/QI/QR/QS/QG`,
`TS`, `MCOSZ`, `FRLAND`, `FROCEAN`). SATSIM has no internal state, so only
import and export are checkpointed.

Assumptions (`radiation.yaml`, `satsim.yaml`, `data-for-satsim.yaml`,
`history.yaml`):

- `USE_SATSIM: 0`, `USE_SATSIM_ISCCP: 1`, all other simulators off, set in
  `radiation.yaml`; `SATSIM_DT: 1200` (= clock step) in `satsim.yaml`. `SATSIM_NCOLUMNS` and
  `SATSIM_POINTS_PER_ITERATION` are left at their defaults (30 at C12, and
  -999 = all points at once).
- Units and standard names of the imports are SatSim's own (as in MAPL2's
  `Satsim_StateSpecs.rc`), so they match the baseline's import metadata.
- In the full model `MCOSZ` comes from SOLAR and the rest come down from
  RADIATION; here `DataForSatsim` supplies all of them.
- No `SATORB` input: the orbit-masked exports are removed in the MAPL3 port.
- Exports are activated through a History collection (`geosgcm_satsim`,
  56 fields) whose start time is after the cap's stop time, so it never
  writes a file. It requests only exports that MAPL2 fills with ISCCP
  alone. The lidar, MISR and MODIS exports would be computed from
  never-written work arrays, so they are not requested.
- The ungridded ISCCP collection (`CLISCCP1-7`, `ISCCP1-7`, size-7
  ungridded dim) is left out: MAPL3 History cannot yet connect to exports
  with ungridded dims (TODO in `history.yaml`).
- SatSim's exports keep MAPL2's undef value (`1e15`).

### `radiation`

The RADIATION container, which creates IRRAD, SOLAR and SATSIM as children
(same `irrad.yaml`, `solar.yaml` and `satsim.yaml` as the standalone cases).
`DataForRadiation` supplies RADIATION's own `PLEINST` plus the union of
IRRAD's, SOLAR's and SATSIM's imports, except `MCOSZ`, which RADIATION
connects from SOLAR to SATSIM. Each is connected once to RADIATION and reaches every child that
imports it; `AERO` comes from FakeGOCART. IRRAD and SOLAR internals are
restarted from `checkpoints/2000-04-14T22:00:00/{IRRAD,SOLAR}_internal.nc`.

Assumptions (`radiation.yaml`, `data-for-radiation.yaml`):

- SatSim runs as in `satsim-sa`: `USE_SATSIM_ISCCP: 1` in `radiation.yaml`
  makes RADIATION create SATSIM and is inherited by it, with `satsim.yaml`
  and `history.yaml` copied from `satsim-sa` (ISCCP only, `SATSIM_DT: 1200`, exports activated
  by the never-writing `geosgcm_satsim` collection). SatSim's exports
  propagate up as `SATSIM/<name>`, so `history.yaml`'s `SATSIM.<name>`
  sources resolve unchanged. `USE_OCEANOBIOGEOCHEM: 0`.
- `DataForRadiation` declares RADIATION's shared imports with the IRRAD/SOLAR
  standard names and units, so SATSIM's imports log standard-name warnings,
  and its `QV`, `QL`, `QI` (units `1`) get a `kg kg-1` -> `1` units
  conversion. `ZLE`, `FRLAND`, `FROCEAN` use SatSim's own metadata.
- `RAD_CONDENSATE_INHOMOGENEITY: 1`, read in RADIATION's `Initialize`.
- IRRAD and SOLAR declare `FCLD` with different standard names; the
  provider uses SOLAR's (`cloud_area_fraction`), so IRRAD's `FCLD` import
  logs a standard-name warning.
- MAPL2's RADIATION import state holds every child import propagated up;
  MAPL3's holds only RADIATION's own (`PLEINST` and the IRRAD/SOLAR exports
  it imports). The children's imports are therefore checkpointed and
  compared through IRRAD's, SOLAR's and SATSIM's own import states.

## Baselines

There are two kinds of baseline: a MAPL2 reference capture, and the MAPL3
checkpoints that ctest compares against.

### MAPL2 reference capture

- **Model**: GEOSgcm **v11.11.0** (MAPL2), with GEOSradiation_GridComp
  **v2.1.0 + fix 1 of
  [#110](https://github.com/GEOS-ESM/GEOSradiation_GridComp/issues/110)**,
  C12L181, run on discover
  (`/discover/nobackup/pchakrab/runs/agcm/capture-radiation-state/C12L181`).
  - #110: SOLAR's `UPDATE_EXPORT` caches the RRTMGP band limits
    (`SOLAR_bands_wavenum`) behind an implicitly saved flag, but the array
    itself is not saved. After the first call, the `TBRBbbRG` (and OBIO)
    exports read uninitialized band limits; a Debug build traps there.
  - Fix 1 deletes the two `have_rrtmgp_wavenums = .true.` lines, so the band
    limits are reloaded on every call. The MAPL3 port has the same change
    (commit `5639a0b`, `TODO` comment in `GEOS_SolarGridComp.F90`), so both
    sides compute these exports the same way.
- **What was captured**: the import, internal and export states of IRRAD,
  SOLAR, SATSIM and RADIATION immediately before and after run phase 1, at
  `2000-04-14T22:00:00`. The files were converted to the new cubed-sphere
  file layout with `GEOS_Restart_Utilities/convertOldCStoNewCS.py`.
- **Run configuration** chosen to match the standalone cases:
  - `AGCM.rc`: `AERO_PROVIDER: none` (default `GOCART2G`).
    `GEOS_ChemGridComp` then provides an empty `AERO` state with
    `implements_aerosol_optics_method = .false.`, which FakeGOCART
    reproduces. A capture made with GOCART2G does not match: SOLAR's
    `FSWN/FSCN/FSWUN/FSCUN/DR*/DF*/FSWBANDN` differ by O(1e-5 to 7e-3) on day
    cells.
  - Radiation resources as in `solar.yaml` / `irrad.yaml` (RRTMGP LW and SW,
    20-minute `SOLAR_DT`/`IRRAD_DT`), and SatSim ISCCP-only with the
    `satsim.yaml` settings.
  - `HISTORY.rc` requesting SOLAR's `[IO]SRB{08..11}RG` and the SATSIM
    ISCCP exports (`geosgcm_satsim`, plus `geosgcm_satsim_ungrd` for the
    ungridded fields). MAPL2 allocates exports only when they are requested.
- **Used for**:
  - The restart files in each case's `checkpoints/2000-04-14T22:00:00/`
    (`DataFor*_internal.nc` from the "before" import states, and
    `IRRAD_internal.nc` / `SOLAR_internal.nc` from the "before" internal
    states).
  - The MAPL2-vs-MAPL3 check done by each case's `cmpchk.py` against the
    "after" states (see below).

### ctest baselines

- `checkpoints/last/*.nc` in each case's regression-data directory: output
  of the MAPL3 port itself (this repo, branch
  `feature/pchakrab/port-satsim-to-mapl3`), generated 2026-10-23 by running
  the cases from the restarts above.
- `run_case.cmake` compares every `*.nc` in the baseline's
  `checkpoints/last` with the run's, using
  `nccmp -dmfgsB --nans-are-equal`, excluding `lons`, `lats`,
  `corner_lons`, `corner_lats`. The comparison must be exact (no
  tolerance).
- The comparison runs **only for IntelLLVM (ifx) builds**. With other
  compilers the case runs but is not compared.

## Regression data

`esma_add_regression_tests` sets the data directory at configure time to
`$LOCAL_REGRESSION_DATA_DIR/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSradiation_GridComp`.
If `LOCAL_REGRESSION_DATA_DIR` is unset, it uses the same path under
`<build>/regression-data/`, synced from S3 by the `irrad/sync_data`
fixture. Layout:

```
GEOSradiation_GridComp/
  ExtData/g5gcm/...                      # RRTMGP v1.8 tables, NRLSSI2.txt (linked as ExtData/)
  <case>/checkpoints/
    2000-04-14T22:00:00/*_internal.nc    # restarts, copied into the run directory
    last/*.nc                            # ctest baseline
```

If a case's data directory is missing, `run_case.cmake` prints
"Regression data not found ... skipping" and the test **passes without
running**.

## Running

```
ctest -L REGRESSION                # all regression tests in the build
ctest -R 'irrad-sa|solar-sa|satsim-sa|radiation'
```

### Comparing against MAPL2 with `cmpchk.py`

After a run, `cmpchk.py` from the run directory compares
`checkpoints/last/<COMP>_<state>.nc` with the MAPL2 "after" capture
(`<COMP>_<state>_after_runPhase1.nc`). It prints, per variable, the
min/max/L2 norm of the absolute difference, and reports separately any
cells that are undefined on one side only. A cell counts as undefined if
it is NaN, `|huge|` (MAPL2's undefined value is `+huge`, MAPL3's
`MAPL_UNDEF` is `-huge`) or either file's `_FillValue`. Grid and dimension variables
are skipped, as are variables missing from either file.

`baseline_dir` is hard-coded to the capture's location on the developer's
machine (`/home/pchakrab/input/radiation/C12L181/after/new-style`); edit it
to point at your copy. The script needs `numpy` and `netCDF4` (e.g. after
`module load ifx-stack`).
