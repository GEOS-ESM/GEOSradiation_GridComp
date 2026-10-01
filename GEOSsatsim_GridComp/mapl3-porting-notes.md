# GEOS_SatsimGridComp MAPL3 porting plan

Status: not started (2026-09-30). `GEOSsatsim_GridComp` is commented
out in the parent container's `alldirs` (`../CMakeLists.txt`) and its
`use` line in `../GEOS_RadiationGridComp.F90` is commented out with
"NOT ported to MAPL3 yet". Nothing in this directory has been touched.

This plan follows the IRRAD/Solar ports (see
`../GEOSirrad_GridComp/mapl3-porting-notes.md` for the MAPL3 API
discoveries/gotchas and `../GEOSsolar_GridComp/mapl3-porting-notes.md`
for the step-by-step template). Only SatSim-specific differences are
spelled out below; everything else is "do what Solar did".

## Reference files
- Template component: `../GEOSsolar_GridComp/GEOS_SolarGridComp.F90`
- Template state spec: `../GEOSsolar_GridComp/Solar_StateSpecs.rc`
- Template build file: `../GEOSsolar_GridComp/CMakeLists.txt`
- Template standalone regression: `../regression/solar-sa/`,
  `../regression/irrad-sa/`, `../regression/cases.txt`
- Parent container: `../GEOS_RadiationGridComp.F90`
- `SATORB` bundle provider (MAPL2): `GEOS_AgcmGridComp.F90` ~L1143
  (`SRC_NAME='SATORB'`/`DST_NAME='SATORB'` connectivity)

## Survey of the MAPL2 code (2026-09-30)
`GEOS_SatsimGridComp.F90` is 5026 lines, all MAPL2:
- `SetServices` (L69-2457): 20 `MAPL_AddImportSpec` (T, QV, FCLD,
  RL, RI, RR, RS, RG, QL, QI, QR, QS, QG, PLE, ZLE, MCOSZ, FRLAND,
  FROCEAN, TS, plus the `SATORB` bundle), 224 `MAPL_AddExportSpec`
  (32 with `UNGRIDDED_DIMS=(/7/)` for the ISCCP/MISR/MODIS 7-bin
  histograms), no internal specs, no FRIENDLYTO/PRECISION/RESTART
  attributes. Private `SatSim_State`/`SatSim_Wrap` set with
  `ESMF_UserCompSetInternalState`. Timers `DRIVER`, `-MISC`, `-COSP`.
  Optional `SatSim.rc` `Masked_Exports::` table parsed with
  `ESMF_Config*` (L2364-2404), then for each row with a `newvar`:
  `MAPL_GetObjectFromGC` -> `MAPL_Get(STATE, ExportSpec=)` ->
  `MAPL_VarSpecGetIndex` -> `MAPL_VarSpecGet(LONG_NAME, UNITS, Dims,
  VLocation, ungridded_dims)` -> dynamic `MAPL_AddExportSpec` cloning
  the spec under `newvar_name` + `MAPL_DoNotDeferExport` (L2405-2445).
  Ends with `MAPL_GenericSetServices`. No Initialize is registered.
- `RUN` (L2465-5022): `MAPL_GetObjectFromGC`, `MAPL_Get(STATE, IM, JM,
  LM, CF, LONS, LATS, INTERNAL_ESMF_STATE, RUNALARM=ALARM)`, returns
  early unless the alarm rings, timer `DRIVER`, calls contained
  `SIM_DRIVER(IM,JM,LM,RC)`. `LONS`/`LATS` are fetched and never used.
- `SIM_DRIVER` (L2799-5020): resources `USE_SATSIM`,
  `USE_SATSIM_ISCCP/MODIS/RADAR/LIDAR/MISR` (integer, default 0),
  `AGCM_GRIDNAME:` (string-parsed for `imsize`/dateline to compute
  `default_Ncolumns = MIN(30, MAX(1, INT(4*5760*4/imsize)))`),
  `SATSIM_NCOLUMNS:`, `SATSIM_POINTS_PER_ITERATION:` (default -999);
  ~250 `MAPL_GetPointer` for imports/exports with ~250 matching pointer
  declarations; `getvistau`/`getirtau` from `gettau`
  (`GEOS_RadiationShared`); COSP calls under timer `-COSP`; masking
  pass at the end (L4884-4990): `ESMF_UserCompGetInternalState`, if
  `nmask_vars > 0`: `ESMF_StateGet(IMPORT, 'SATORB', BUNDLE)`,
  `MAPL_Get(STATE, ExportSpec=)`, `MAPL_VarSpecGetIndex/Get(Dims,
  ungridded_dims)`, `ESMFL_BundleGetPointerToData(BUNDLE, mask_name,
  ptr_mask)`, sets export to `MAPL_UNDEF` where `ptr_mask ==
  MAPL_UNDEF` (2D, 2D+ungridded and 3D cases). Note the `nullify`
  inside the ungridded loop at L4931-4932 looks like an existing bug.
- Macro counts: `VERIFY_` 536, `RETURN_` 5, `__RC__`/`__STAT__` used,
  `MAPL_UNDEF` 74, `MAPL_Am_I_Root` 4, `DEBUG_GC=.FALSE.` debug prints.
- The COSP/quickbeam/icarus/MISR/MODIS numerical sources
  (`cosp*.F90`, `modis_simulator.F90`, `icarus.f`, `isccp_cloud_types.f`,
  `MISR_simulator.f`, `scops.f`, `geos5_*.f90`, `actsim/`, `llnl/`,
  `quickbeam/`, `cmor/`) are MAPL-free and need no porting.
- No `SatSim.rc` exists anywhere in the tree; the masked-exports path
  has probably never been exercised in this repo.

## Steps

1. **Reformat the GridComp skeleton first.** DONE 2026-10-01.
   - DONE - `codee format --verbose ~/config/codee-format
     GEOS_SatsimGridComp.F90` ran clean, no `!CODEE_MASK!` needed
     (5026 -> 5122 lines, same as the 2026-09-30 trial).
   - DONE - macro conversion (regex script, then hand edits): 530
     `RC=STATUS)`+`VERIFY_` pairs -> `_RC)`, 16 `stat=STATUS`/`__STAT__`
     -> `_STAT`, 4 `RETURN_(ESMF_SUCCESS)` -> `_RETURN(_SUCCESS)`, the 2
     positional-status `ESMF_UserComp{Set,Get}InternalState` calls ->
     `_VERIFY(status)`, `STATUS` -> `status`. Fixed a pre-existing
     missing `VERIFY_` after the `PARASOLREFL0` export spec. Kept the two
     intentional soft-fail `if (status == ESMF_SUCCESS)` checks in the
     `Masked_Exports::` parsing (optional table / optional third column).
     Dropped `IAm`/`COMP_NAME` declarations + the three
     `ESMF_GridCompGet(NAME=COMP_NAME)` traceback lookups, all
     `!BOP/!EOP/!BOS/!EOS` blocks, and replaced the Moist copy-paste
     module description with a one-paragraph SatSim one. File is now
     4513 lines. Structural counts (subroutines, `#if/#endif`, 244
     AddSpec, 250 GetPointer) match the original.

2. **`Satsim_StateSpecs.rc`** DONE 2026-10-01.
   - DONE - `Satsim_StateSpecs.rc` written (generated from the MAPL2
     spec calls by a throwaway script, then validated with
     `MAPL_GridCompSpecs_ACG.py`): 19 field imports (`PLE`/`ZLE` are
     `VLOC=E`; the `SATORB` bundle import is excluded, ACG ITEMTYPE
     only supports F/V) and 222 static exports (30 with `UNGRIDDED:
     ungrd_7`, `PARASOLREFL0` with `ungrd_parasol_nrefl`). The earlier
     "20/224" counts included `SATORB` and the 2 dynamic masked-export
     specs. Imports carry an `ALIAS` column so the generated pointers
     keep the historical `SIM_DRIVER` names (`RL->RDFL`, `QL->QLTOT`,
     `QR->QRTOT`, ...). ACG output: 19 + 222 `MAPL_GridCompAddSpec`,
     241 `MAPL_StateGetPointer` = the 241 removed `MAPL_GetPointer`
     calls; pointer names and ranks cross-checked against the removed
     hand declarations (no mismatches).
   - DONE - `CMakeLists.txt`: `mapl_acg(${this} Satsim_StateSpecs.rc
     IMPORT_SPECS EXPORT_SPECS GET_POINTERS DECLARE_POINTERS)`.
   - DONE - `SetServices`: declared `type(MAPL_UngriddedDim) :: ungrd_7,
     ungrd_parasol_nrefl`, set them (`MAPL_UngriddedDim(7,
     name='bins7', units='1')`, `MAPL_UngriddedDim(PARASOL_NREFL,
     name='parasol_nrefl', units='1')`) and replaced the 19+222 spec
     calls with `#include "Satsim_Import___.h"` /
     `"Satsim_Export___.h"`. `SATORB` kept as a hand-written MAPL2
     `MAPL_AddImportSpec` for step 3; the 2 dynamic masked-export specs
     untouched for step 4.
   - DONE - `Run`: all export pointer declarations replaced by
     `#include "Satsim_DeclarePointer___.h"` (which also declares the
     imports; `SIM_DRIVER` host-associates them, its duplicate import
     declarations were removed). `SIM_DRIVER`: the 241 `MAPL_GetPointer`
     calls replaced by `#include "Satsim_GetPointer___.h"`. Only the
     dynamic masked-export `MAPL_GetPointer` calls remain (step 4).
     File is now 2112 lines (from 4513).

3. **`SATORB` import bundle** - the only non-field import
   (`ESMF_StateGet(IMPORT, 'SATORB', BUNDLE)`). Decide how a bundle
   import is declared in MAPL3 (StateSpecs or a manual spec) and how
   the parent connects it; `ESMFL_BundleGetPointerToData` needs a
   MAPL3 equivalent (`ESMF_FieldBundleGet` + `ESMF_FieldGet(farrayPtr)`
   or a MAPL3 bundle-pointer helper).

4. **Masked exports (`SatSim.rc` `Masked_Exports::` table)** - the
   hard part. Read the table from the component HConfig (`satsim.yaml`,
   a sequence of `{export, mask, newvar}` maps; HConfig API in IRRAD
   notes) and add the new exports via the MAPL3 dynamic-spec API (as
   IRRAD did for its runtime-named RATS_DIAGNOSTICS exports - dynamic
   fields cannot go in the YAML). Copying dims/units/ungridded from an
   existing spec needs a MAPL3 way to introspect a registered spec, or
   re-declare the copy's metadata from a small lookup of the handful of
   shapes SatSim has (2D, 2D+7, 3D). `MAPL_DoNotDeferExport` has no
   MAPL3 counterpart - check whether it is still needed.

5. **Private state** - `SatSim_State`/`SatSim_Wrap` +
   `ESMF_UserCompSet/GetInternalState` -> `_SET_NAMED_PRIVATE_STATE` /
   `_GET_NAMED_PRIVATE_STATE` (setter and both getters together). The
   state only holds the mask table, so it could instead be re-read from
   HConfig in `Run` and the private state dropped.

6. **Resources** - 9 `MAPL_GetResource(MAPL, ..., LABEL="X:")` ->
   `MAPL_GridCompGetResource(gc, "X", ..., default=..., _RC)`. Derive
   `imsize` for the `SATSIM_NCOLUMNS` default from the grid
   (`MAPL_GridCompGet` / `ESMF_GridGet` global dims) instead of
   string-parsing `AGCM_GRIDNAME`. `USE_SATSIM*` integer flags: keep as
   integers with `default=0`, or convert to logicals.

7. **`Run` plumbing** - drop `MAPL_MetaComp` / `MAPL_GetObjectFromGC` /
   `MAPL_Get`: `IM/JM/LM` via `MAPL_GridCompGet`/grid; delete the unused
   `LONS/LATS`; `RUNALARM` -> create a `satsim_alarm` in a real
   `Initialize` using Solar's `RUN_AT_INTERVAL_START` /
   `REFERENCE_DATE` / `REFERENCE_TIME` block with `SATSIM_DT`, plus
   `_ASSERT` `DT /= 0` and `mod(DT, clock dt) == 0` (IRRAD notes
   2026-09-29). `MAPL_TimerAdd/On/Off` -> `MAPL_GridCompTimerStart/Stop
   (gc, name, _RC)`. Remove `DEBUG_GC` and the `write(*,*)` debug prints
   (or route through the logger).

8. **Edge remaps** - `PLE`/`ZLE` are `VLOC=E`; MAPL3 pointers are
   1-based but `SIM_DRIVER` indexes `PLE(:,:,0:LM)` and declares
   `PLE2D(IM*JM,0:LM)` locals. Apply the contiguous 0-based rank-remap
   trick from Solar step 11 / IRRAD notes, or shift the indexing.

9. **Entry-point signatures** - `SetServices(gc, rc)` /
   `Initialize(gc, import, export, clock, rc)` / `Run(gc, import,
   export, clock, rc)`, `MAPL_GridCompSetEntryPoint` for both phases,
   MAPL3-style `SetServices` ending (no `MAPL_GenericSetServices`). The
   container's `use ... only: satsimSetServices => SetServices` alias
   should suffice, as for IRRAD/SOLAR.

10. **Build wiring** - `GEOSsatsim_GridComp/CMakeLists.txt`: fix
    `DEPENDENCIES` (needs `GEOS_RadiationShared` for `gettau`, `MAPL`,
    `GEOS_Shared`), `TYPE SHARED` to match siblings. In the parent:
    uncomment `GEOSsatsim_GridComp` in `alldirs` and the `use
    GEOS_SatsimGridCompMod` line; `MAPL_GridCompAddChild(gc, "SATSIM",
    satsimSetServices, "satsim.yaml", _RC)` guarded by `USE_SATSIM` as
    in MAPL2; connect the imports (most already flow through Radiation
    for IRRAD/SOLAR; `SATORB`, `FROCEAN`, `TS`, `MCOSZ`, hydrometeor
    `Q*`/`R*` need checking) and re-export the 224+ exports via
    `MAPL_GridCompReexport` - loop over the ACG-generated list rather
    than 224 hand lines; check whether MAPL3 has a "reexport all of
    child" helper. Add `satsim.yaml` (`SATSIM_DT`, `USE_SATSIM_*`,
    `SATSIM_NCOLUMNS`, optional `Masked_Exports`).

11. **Compile with ifx** and fix fallout. Only the GridComp should
    break; the fixed-form COSP `.f` sources are untouched.

12. **Standalone regression `regression/satsim-sa`** - clone
    `irrad-sa` (`cap_*.yaml`, `mapl.yaml`, `root.yaml`, `data-*.yaml`
    provider, `cmpchk.py`). The fake provider must supply all 20 imports
    plus a `SATORB` bundle with the mask fields, using dictionary
    `standard_name`s (Solar notes RR lesson). Generate a MAPL2 baseline
    (SatSim has never been run-tested in this repo, same as IRRAD was)
    and add the case to `regression/cases.txt` for ctest.

13. **Style cleanup afterwards** (low priority, like Solar step 16):
    lowercase keywords, collapse column-aligned declarations, remove the
    commented-out `frac_outinv` / `congvec` / `load_mie_table` cruft,
    and fix the `nullify`-inside-loop bug in the masking pass (L4931)
    once a test covers it.

Effort concentrates in steps 2 (mechanical but ~2300 lines), 4 (needs a
MAPL3 design decision for dynamic export cloning) and 12 (test data and
baseline).
