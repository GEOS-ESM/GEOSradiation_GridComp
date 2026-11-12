# GEOSradiation_GridComp: guide for agents

**Describes:** the MAPL3 line, branch `release/MAPL-v3`; verified against `GEOSradiation_GridComp` @ `856688c` (2026-11-10).
**Placement:** `src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSradiation_GridComp` in the GEOSgcm fixture (`components.yaml`); other fixtures may place it elsewhere. Paths below are relative to this repository.

> **Built and ported on this line.** `CMakeLists.txt` builds the `RADIATION` composite and `GEOS_RadiationShared`, `GEOSsolar_GridComp`, `GEOSirrad_GridComp` and `GEOSsatsim_GridComp`, all ported to MAPL3, and adds `regression/`. `GEOS_RadiationShared/@RRTMGP` is a separate nested repository; edit it only there. The Physics library that would add `RADIATION` (`GEOS_PhysicsGridComp.F90` in `GEOSgcm_GridComp`) is not compiled on this line, so `RADIATION` runs only in the regression cases; the comment "restricted to IRRAD only" on this repository in Physics' `CMakeLists.txt` predates the ports.

## Rules for agents

- Read [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) (what Radiation is and how it is put together: structure, refresh and every-step split, scheme selection, contracts, invariants, the background to the traps below) before changing children, connections, exports, schemes or configuration. Update it when the structure, a contract or an invariant changes; it describes only what changes rarely, so keep status and to-dos out of it.
- Note bugs and oddities; do not fix them unless asked. Record a bug or an open question as a GitHub issue, not in a document; a document may link to the issue.
- Do not edit the generated `*___.h` files; edit the `*_StateSpecs.rc` table (`Radiation_StateSpecs.rc`, `GEOSirrad_GridComp/Irrad_StateSpecs.rc`, `GEOSsolar_GridComp/Solar_StateSpecs.rc`, `GEOSsatsim_GridComp/Satsim_StateSpecs.rc`). Only state the table cannot express is added in `SetServices`.
- When porting code to MAPL3, follow "Porting to MAPL3" in the `GEOSgcm_GridComp` `AGENTS.md` (`../../../AGENTS.md` in the fixture).
- Keep this file current: update it when status, a rule or a trap changes, and update its "verified against" commit.
- Cite files and routines, never line numbers. Label unchecked statements *Inferred:* or *Unverified:*.
- Error handling in code you write or touch uses the MAPL3 macros (`@MAPL/include/MAPL_ErrLog.h`): `_RC` on calls, `_STAT` on `allocate`, `_ASSERT(cond, msg)` with a message, and `_RETURN(_SUCCESS)` at every successful exit. Convert `RC=STATUS`/`VERIFY_`/`ASSERT_`/`RETURN_` in any routine you change. `_RC` expands to two statements (`rc=status);_VERIFY(status`), so it cannot end a one-line `if`; use a block `if`. Do not add decorative `!----` separator lines.

## Status and to-dos

- `IRRAD`'s `SetServices` does not read the scheme keys: it sets `USE_RRTMGP = .true.` ("TODO: determine from configuration; default to RRTMGP for now"), so its scheme-dependent specs are declared for RRTMGP whatever the keys say. `Run` reads the keys.
- `IRRAD` and `SOLAR` call `set_inhomogeneity` in their own `Initialize`, for standalone runs ("TODO: eventually remove this"); the composite already sets it.
- `SOLAR` keeps fix 1 of GEOSradiation_GridComp#110 (band limits reloaded on every call), with a `TODO` in `UPDATE_EXPORT`.
- `SOLARAvrg` and `IRRADAvrg` (MAPL2 import averaging) are not read on this line; a refresh sees the imports of the step on which its alarm rings.
- `SATSIM`: the MAPL2 `SATORB` orbit masking is removed, and the ungridded ISCCP exports cannot yet be written by MAPL3 History.
- Porting history: `mapl3-porting-notes.md` at the root and in `GEOSsolar_GridComp/`, `GEOSirrad_GridComp/`, `GEOSsatsim_GridComp/`.

## Known traps

- **A new scheme touches several places:** both children (`SetServices` of `SOLAR`, `Run` of both), `rad_utils.F90`, the band counts and the `NUM_BANDS` check.
- **`NUM_BANDS` is not derived.** Changing a scheme key without changing `NUM_BANDS` stops both children in their first `Run` ("Total number of radiation bands is inconsistent!"); the aerosol provider's optics tables must have the matching bands.
- **A change to clouds, humidity or albedo reaches the heating only at the next refresh**; between refreshes only insolation and surface temperature act.
- **The composite cannot read a child export directly.** Using one needs an import in `Radiation_StateSpecs.rc` and a `MAPL_GridCompAddConnection` to `"<self>"` in `SetServices`; otherwise the child leaves it unfilled.
- **Re-bound Edge pointers before using them** with the `0:LM` convention (as `Run`, `Update_Flx` and the child `Run` routines do); MAPL3 gives them `1:LM+1`.
- **Two code paths per flux.** A new flux export needs a stored internal and a line in the every-step update (`UPDATE_EXPORT` in `SOLAR`, `Update_Flx` in `IRRAD`) as well as in the refresh.
- **New work for an export must be guarded like the rest:** done only when the export is requested.
- **`SOLAR` packs its imports and internals for load balancing** (`SORADCORE`): every import is packed by name, only 2-D horizontal and 3-D imports are accepted (`AERO` is the one non-field import allowed), internals are packed with name-based cases (`INT_VARS_1`), and night points take the internal's `FILL` value, so keep the `FILL` column in step with the code.
- **`IRRAD` sign conventions:** upward fluxes are stored negated (net downward positive); `DFDTS` is the derivative of the negated upward flux.
- **Band-resolved outputs exist only for some bands** (`band_output_supported` in `SOLAR` and `IRRAD`); adding a band needs the table, the internals and the exports.
- **`SOLAR_DT` is not checked:** `IRRAD` and `SATSIM` assert that their interval is a positive multiple of the time step; `SOLAR` does not.
- **`CALLED_LAST: 0` is inconsistent with `SOLAR`'s alarm**, which a source comment says is always set for the last-step case (`Run`, "Determine calling sequence").
- **RRTMGP shortwave needs the NRLSSI2 solar-cycle file** (`SOLAR_CYCLE_FILE_NAME` not `/dev/null`, `USE_NRLSSI2`); `SORADCORE` asserts it.
- **`ENABLE_SOLAR_RADVAL` changes `SOLAR`'s restart contents:** restarts from builds with and without it are not interchangeable.

## Build and test

- Builds only as part of the GEOSgcm fixture (see its `README.md`); this repository has no standalone build.
- Regression tests (label `REGRESSION`): `regression/`, cases `irrad-sa`, `solar-sa`, `satsim-sa`, `radiation` (`cases.txt`). `regression/README.md` describes the cases, their data, what is compared, and their run time. Each case directory has a `logging.yaml` (keep it in new cases). A case whose data directory is missing passes without running, and the comparison is done only for IntelLLVM (ifx) builds (`run_case.cmake`).
- Pull requests for this line target `release/MAPL-v3`, not `develop` (the MAPL2 line), and need a 0-diff label (`.github/workflows/enforce-labels.yml`). CI builds the fixture at `release/MAPL-v3` and runs the regression tests (`.circleci/config.yml`, `.github/workflows/workflow.yml`). Code owners: the radiation team, and the CMake team for `CMakeLists.txt` files (`.github/CODEOWNERS`).
