# AGENTS.md

**Audience:** coding agents and developers working in the `GEOSradiation_GridComp` repository.
**Describes:** the MAPL2 line, branch `develop`.
**Verified against:** `GEOSradiation_GridComp` @ `0b402f6` (2026-11-05).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `0f2e0fb`).
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

## What this is

The radiation part of the GEOS atmospheric physics: the `RADIATION` gridded component
(`GEOS_RadiationGridComp.F90`) and its children, shortwave `SOLAR`, longwave `IRRAD` and
the optional satellite simulator `SATSIM`, plus shared code. In the GEOSgcm fixture it is at
`src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSradiation_GridComp`
(`components.yaml`, key `GEOSradiation_GridComp`) and is the last child run by Physics. The
`GEOS_RadiationShared/@RRTMGP` directory is a separate nested repository (`RRTMGP`); edit it
only in its own repository.

## Documents

Read `docs/ARCHITECTURE.md` first, then the file for the directory you touch.

| Document | Content |
|---|---|
| [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) | Placement, directory map, the composite, scheme selection, refresh versus every-step split, aerosol optics from Chemistry, interfaces to Physics siblings, build |
| [`docs/decisions.md`](docs/decisions.md) | Previous architectural decisions (pattern, likely reason) |
| [`GEOSsolar_GridComp/ARCHITECTURE.md`](GEOSsolar_GridComp/ARCHITECTURE.md) | Shortwave child: Chou-Suarez, RRTMG, RRTMGP, load balancing, solar forcing |
| [`GEOSirrad_GridComp/ARCHITECTURE.md`](GEOSirrad_GridComp/ARCHITECTURE.md) | Longwave child: three schemes, surface-temperature update, RATS diagnostics |
| `README.md`, `CONTRIBUTING.md` | Project basics, contribution rules |

The satellite simulator (`GEOSsatsim_GridComp/`) and the shared code
(`GEOS_RadiationShared/`) are described in sections 8 and 9 of `docs/ARCHITECTURE.md`.
Parent documents (Physics, GEOSgcm_GridComp, fixture) are linked from the header of
`docs/ARCHITECTURE.md`.

## Rules

- Note bugs you see; do not fix them unless asked. Record them under "Observations" in the
  relevant document.
- Keep the documents current: if you change a child list, run order, connection, resource key,
  scheme choice or build option described there, update the document in the same change. Cite
  files and routines, not line numbers.
- Do not edit generated files (`*_Import___.h`, `*_Export___.h`, `*_Internal___.h`, ...); edit
  `GEOSirrad_GridComp/Irrad_StateSpecs.rc`. `SOLAR` and `SATSIM` have no ACG file and declare
  their state in `SetServices`.
- Code owners are in `.github/CODEOWNERS`: the radiation team, and the CMake team for
  `CMakeLists.txt` files.

## Build and test

This repository does not build alone; build the GEOSgcm fixture (see its `README.md`). Build
structure: `docs/ARCHITECTURE.md` section 10. No unit-test target was found. CI builds the
fixture and runs short model runs (`.circleci/config.yml`, `.github/workflows/workflow.yml`);
pull requests need a 0-diff label (`.github/workflows/enforce-labels.yml`).
