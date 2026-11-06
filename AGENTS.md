# AGENTS.md

**Audience:** coding agents and developers working in the `GEOSradiation_GridComp` repository.
**Describes:** the MAPL3 line, branch `release/MAPL-v3`.
**Verified against:** `GEOSradiation_GridComp` @ `13f8eb8` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `f50622d`).
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

## What this is

> **Build status on this line:** built, with standalone regression cases. On
> `release/MAPL-v3`, `CMakeLists.txt` builds the `RADIATION` composite and its
> subdirectories `GEOS_RadiationShared`, `GEOSsolar_GridComp`, `GEOSirrad_GridComp`
> and `GEOSsatsim_GridComp`, and adds `regression/`; the code is the MAPL3 port
> (history in the `mapl3-porting-notes.md` files). The Physics `CMakeLists.txt`
> comment on this repository ("restricted to IRRAD only") predates the ports (*Inferred*); the
> Physics library itself is not compiled on this line, so the MAPL2 Physics code that
> adds `RADIATION` is not part of the build.

The radiation part of the GEOS atmospheric physics: the `RADIATION` gridded component
(`GEOS_RadiationGridComp.F90`) and its children, shortwave `SOLAR`, longwave `IRRAD` and
the optional satellite simulator `SATSIM`, plus shared code. In the GEOSgcm fixture it is at
`src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSradiation_GridComp`
(`components.yaml`, key `GEOSradiation_GridComp`) and is the last child run by Physics (in the MAPL2 Physics code, which this line does not compile). The
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
| [`regression/README.md`](regression/README.md) | The four standalone regression cases (assumptions, data, what is compared) |
| `mapl3-porting-notes.md` (root and in the three child directories with code ports) | Porting history and MAPL3 findings; refer to them, they are not repeated in the documents above |
| `README.md`, `CONTRIBUTING.md` | Project basics, contribution rules |

The satellite simulator (`GEOSsatsim_GridComp/`, no document of its own) and the shared code
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
  the ACG input: `Radiation_StateSpecs.rc`, `GEOSirrad_GridComp/Irrad_StateSpecs.rc`,
  `GEOSsolar_GridComp/Solar_StateSpecs.rc`, `GEOSsatsim_GridComp/Satsim_StateSpecs.rc`. State
  that the ACG file cannot express (the `AERO` state item, the `RATS_DIAGNOSTICS` fields,
  the solar band exports built in loops) is added in `SetServices`.
- Code owners are in `.github/CODEOWNERS`: the radiation team, and the CMake team for
  `CMakeLists.txt` files.

## Build and test

This repository does not build alone; build the GEOSgcm fixture (see its `README.md`). Build
structure: `docs/ARCHITECTURE.md` section 10. The tests are the four standalone regression cases in
`regression/` (`irrad-sa`, `solar-sa`, `satsim-sa`, `radiation`; ctest label `REGRESSION`; they need
their data directory; `regression/README.md`). CI builds the fixture at `release/MAPL-v3` and runs
them (`.circleci/config.yml`, `.github/workflows/workflow.yml`); pull requests need a 0-diff label
(`.github/workflows/enforce-labels.yml`).
