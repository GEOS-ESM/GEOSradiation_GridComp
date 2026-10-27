# GEOS_SatsimGridComp MAPL3 porting plan

Status (2026-10-22): in progress on `feature/pchakrab/port-satsim-to-mapl3`.
Steps 1-11 and 13 are done (step 4 by removing the masked exports).
SatSim is built as a RADIATION child (behind `USE_SATSIM`) and compiles;
it has not been run. Open: step 12 (the `satsim-sa` regression).
Step 14 (SatSim writes into its radius imports) DONE 2026-10-26.

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

3. **`SATORB` import bundle** DONE 2026-10-01. Commented out since
   2026-10-22, and removed 2026-10-26, with the masked exports, its only
   consumer (step 4).
   - `MAPL_AddImportSpec(..., DATATYPE=MAPL_BundleItem)` ->
     `MAPL_GridCompAddSpec(gc, state_intent=ESMF_STATEINTENT_IMPORT,
     short_name='SATORB', standard_name='Satellite_orbits', dims="xy",
     vertical_stagger=MAPL_VERTICAL_STAGGER_NONE, units='days',
     itemtype=MAPL_STATEITEM_FIELDBUNDLE, _RC)`. This mirrors the
     producer's export spec in
     `src/Shared/@MAPL/gridcomps/orbit/MAPL_OrbGridCompMod.F90` exactly;
     the ACG still cannot express bundles, so it stays hand-written.
   - `ESMFL_BundleGetPointerToData(bundle, name, ptr, _RC)` ->
     `MAPL_FieldBundleGetPointer(bundle, trim(name), ptr, _RC)`
     (`infrastructure/field_bundle/FieldBundleGetPointer.F90`, exported
     from that directory's `API.F90`). Same semantics: it nullifies the
     pointer when the field is not `ESMF_FIELDSTATUS_COMPLETE`, so the
     existing `associated()` guards still hold. `trim` is needed because
     `self%mask_name` is `ESMF_MAXSTR`-padded and the lookup is by exact
     field name.
   - `ESMF_StateGet(import, 'SATORB', bundle, _RC)` is unchanged - that
     is plain ESMF and still correct in MAPL3.
   - The ORBIT producer is already fully MAPL3 (private-state macros,
     `MAPL_GridCompSetEntryPoint`, colon-free resources) and is built as
     the `MAPL.orbit` target, so no work is needed on that side.
   - **Deferred to step 10 (connectivity).** In MAPL2 the only explicit
     wiring was `MAPL_AddConnectivity(GC, SRC_NAME='SATORB',
     DST_NAME='SATORB', SRC_ID=ORB, DST_ID=PHYS)` in
     `GEOS_AgcmGridComp.F90`; `SATORB` reached SATSIM by MAPL2's
     automatic propagation of unsatisfied child imports up the tree
     (SATSIM -> RADIATION -> PHYSICS). **CORRECTED 2026-10-22**: MAPL3
     does propagate unconnected child imports up the tree
     (`StateRegistry_Propagation_smod.F90`; relied on by SU/DU/SS's
     `TRADV` imports, see the SU2G porting notes step 4), so no
     per-level `MAPL_GridCompReexport` of `SATORB` should be needed; a
     single connection from ORBIT at the level where both are visible
     suffices. The blocker is that `GEOS_AgcmGridComp.F90` is still
     MAPL2, so there is no live ORBIT instance to connect to yet.
   - Mask semantics, for reference when building the step 12 fake
     provider: ORBIT sets `field = 1.0` inside the satellite swath and
     `MAPL_UNDEF` everywhere else, which is what SatSim's
     `where (ptr_mask == MAPL_UNDEF)` test keys on. Bundle field names
     are the instrument names from the `Nominal_Orbits::` table in
     `MAPL_OrbGridComp.rc` (`MODIS_T`, `CALIPSO`, `CLOUDSAT`, ...), and
     those are exactly the names `SatSim.rc`'s `Masked_Exports::` table
     puts in its mask column.

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
   - State (2026-10-22): still all MAPL2. `SetServices` reads
     `SatSim.rc` via `ESMF_GridCompGet(CONFIG=)` + `ESMF_Config*`, then
     `MAPL_GetObjectFromGC` -> `MAPL_VarSpecGet` -> `MAPL_AddExportSpec`
     + `MAPL_DoNotDeferExport`; the masking pass (end of `SIM_DRIVER`)
     uses `MAPL_VarSpecGet` and `MAPL_GetPointer`.
   - Spec metadata cannot be introspected at `SetServices` time (no
     `MAPL_VarSpec` in MAPL3; Solar step 12). At `Run` time the masking
     pass can read it from the field instead: `MAPL_FieldGet(field,
     ungridded_dims=, num_levels=)`, as Solar's `SolarFieldGetDims`
     does (`ungridded_dims` excludes the vertical; use `num_levels > 0`
     to tell `xyz` from `xy`+ungridded). Remaining
     `MAPL_GetPointer` -> `MAPL_StateGetPointer`.
   - Options considered: (1) `Masked_Exports` list in `satsim.yaml`, with
     the clone's units/dims/ungridded from a small lookup in SatSim;
     (2) the same, with that metadata carried in each YAML entry;
     (3) drop the feature (no `SatSim.rc` exists anywhere, so it has
     never run here) and leave masking to downstream tools.
   - **DONE 2026-10-22 - disabled, not ported.** The masked-export code
     is commented out (not deleted) with a note to port it to MAPL3 before
     re-enabling: the `SatSim.rc` parsing and spec cloning in
     `SetServices`, the masking pass at the end of `SIM_DRIVER`
     (including its `_GET_NAMED_PRIVATE_STATE`), and the locals only they
     use. The `ESMF_GridCompGet(gc, CONFIG=CF)` call was deleted: `CF`
     was recreated by `ESMF_ConfigCreate` before any use. No MAPL2 calls
     remain outside comments. The `SATORB` bundle import (step 3) is
     commented out too, under the same note: its only consumer was the
     masking pass, and keeping it would mean wiring ORBIT (step 10) and
     faking a bundle in `satsim-sa` (step 12) for nothing. Re-enable it
     with the masked exports. Still present but unused: the
     `SatSim_State` private state (only the mask table,
     `nmask_vars = 0`). Options (1)/(2) above remain the way to port it.
   - **REMOVED 2026-10-26.** The commented-out masked-export code (the
     `SATORB` spec, `SatSim.rc` parsing and spec cloning in `SetServices`,
     the masking pass in `SIM_DRIVER`, and their commented locals) is
     deleted, along with the `SatSim_State` private state (type,
     `PRIVATE_STATE`, `_SET/_GET_NAMED_PRIVATE_STATE` in `SetServices`),
     which held only the mask table. A note in `SetServices` points to
     commit `0a0936c` (last version with the code, commented out) and to
     `release/MAPL-v3` at `7862c74` (the live MAPL2 version).
   - `MAPL_UNDEF` changed value: MAPL2 `1.0e15`
     (`base/Base/Base_Base.F90:78`), MAPL3 `-huge(1.)`
     (`utils/Constants/InternalConstants.F90:18-19`). (MAPL2's `+huge`
     was `MAPL_UNDEFINED_REAL`, a different constant.) The masking test
     `ptr_mask == MAPL_UNDEF` is commented out with the masked exports.
     The radius defaults in `SIM_DRIVER` (~L540-557, `where (RDFL ==
     ...) RDFL = 14.e-6` and the same for `RDFI`, `RDFR`, `RDFS`,
     `RDFG`) key on the undef value of the radius imports, which come
     from MAPL2-era data holding `1e15` (Solar step 18). With
     `MAPL_UNDEF` they never matched. DONE 2026-10-22: they now test
     against a module parameter `MAPL2_UNDEF = 1.0e15`, as IRRAD does
     (`GEOS_IrradGridComp.F90:159`, ~L1021-1025). Every other live
     `MAPL_UNDEF` (the `#ifdef USE_MAPL_UNDEF` blocks that turn COSP's
     `R_UNDEF` and out-of-range values into undef, ~L882-983, and the
     radar mean, ~L1419-1430) is now `MAPL2_UNDEF` too, so SatSim's
     exports carry `1e15` as in MAPL2, like IRRAD's. The radar test at
     ~L1419 depends on ~L970 writing the same constant; both changed
     together. The commented-out masked-export code still says
     `MAPL_UNDEF` (ORBIT's mask would be MAPL3).

5. **Private state** DONE 2026-10-01.
   - Kept the private state (the mask table is parsed once in
     `SetServices`; re-reading HConfig in `Run` would duplicate the
     parse for no gain) but dropped the hand-rolled `SatSim_Wrap`
     type - `_SET/_GET_NAMED_PRIVATE_STATE` declare their own wrapper.
   - Added `character(*), parameter :: PRIVATE_STATE = "SatSim_state"`
     at module scope (Solar's `PRIVATE_STATE` pattern).
   - `SetServices`: `allocate(self)` + `wrap%PTR => self` ->
     `_SET_NAMED_PRIVATE_STATE(GC, SatSim_State, PRIVATE_STATE)`
     immediately followed by `_GET_NAMED_PRIVATE_STATE(..., self)`;
     the trailing `ESMF_UserCompSetInternalState` is gone (the macro
     allocates and attaches in one shot, so the set has to move to the
     top where `self` is first used).
   - `SIM_DRIVER`: `ESMF_UserCompGetInternalState` + `self => wrap%PTR`
     -> `_GET_NAMED_PRIVATE_STATE(GC, SatSim_State, PRIVATE_STATE,
     self)`.
   - Gave `nmask_vars` a `= 0` default: the macro-allocated state is
     no longer zeroed by the old `else self%nmask_vars = 0` path only.
   - The macros come from `MAPL_private_state.h`, which `MAPL.h`
     includes; no new include needed. (The file now includes `MAPL.h`
     directly; `MAPL_Generic.h` is just `#include "MAPL.h"`, step 13.)

6. **Resources** DONE 2026-10-01.
   - All 9 `MAPL_GetResource(MAPL, x, LABEL="X:")` ->
     `MAPL_GridCompGetResource(GC, "X", x, default=..., _RC)` (note:
     no trailing colon on the MAPL3 key). Kept `USE_SATSIM*` as
     integers with `default=0` - they are summed
     (`use_satsim + use_satsim_isccp > 0`) in ~15 places, so converting
     to logicals would mean rewriting all of those.
   - Dropped `AGCM_GRIDNAME` entirely. `imsize` (only used to pick the
     `SATSIM_NCOLUMNS` default) now comes from the grid:
     `MAPL_GridCompGet(GC, grid=esmfgrid)` ->
     `MAPL_GridGetGlobalCellCountPerDim`, then
     `imsize = 4*IM_World` if `JM_World == 6*IM_World` else `IM_World`.
     That is the same cubed-sphere rule MAPL3's own Orbit GridComp uses
     (`MAPL_OrbGridCompMod.F90` L338-345) and reproduces the MAPL2
     `dateline == 'CF'` branch without string parsing. Removed the
     `GRIDNAME`/`imchar`/`dateline`/`nn` locals.
   - Removed the `MAPL_GetObjectFromGC(GC, MAPL, ...)` call and the
     `type(MAPL_MetaComp), pointer :: MAPL` local from `SIM_DRIVER`;
     the two in `SetServices`/`Run` stay until steps 4 and 7.

7. **`Run` plumbing** DONE 2026-10-01.
   - Dropped the `MAPL_MetaComp`/`MAPL_GetObjectFromGC`/`MAPL_Get` block
     from `Run`. `LM` now comes from
     `MAPL_GridCompGet(gc, grid=esmfgrid, num_levels=LM, _RC)` and
     `IM/JM` from `MAPL_GridGet(esmfgrid, IM=IM, JM=JM, _RC)`, matching
     IRRAD's `Run`.
   - Deleted the unused `LONS`, `LATS`, `CF` and `internal` locals. The
     MAPL2 code requested `INTERNAL_ESMF_STATE` but never used it;
     SatSim has no internal state spec.
   - `esmfgrid` is declared once in `Run` and host-associated into
     `SIM_DRIVER`, so the step 6 `MAPL_GridCompGet` call there (for the
     `imsize` default-`Ncolumns` calculation) was removed as redundant.
   - The generic `RUNALARM` is replaced by
     `ESMF_ClockGetAlarm(clock, alarmname="satsim_alarm", alarm=satsim_alarm, _RC)`,
     looking up the alarm created in `Initialize` (step 9). The early
     `_RETURN(_SUCCESS)` when it is not ringing is unchanged.
   - `MAPL_TimerAdd`/`MAPL_TimerOn`/`MAPL_TimerOff` -> `MAPL_GridCompTimerStart`/
     `MAPL_GridCompTimerStop(gc, name, _RC)`. MAPL3 timers do not need
     to be registered, so the three `MAPL_TimerAdd` calls in
     `SetServices` are gone. Note the `-MISC` timer was only ever added,
     never started/stopped, so it disappears entirely.
   - Removed `DEBUG_GC` and the `MAPL_AM_I_Root() .and. DEBUG_GC`
     `write(*,*)` block around the `COSP` call. The flag was
     hard-wired to `.false.`, so the prints were dead code.
   - Carry-over: the masked-export loop at the end of `SIM_DRIVER` still
     needs the MAPL2 generic state to introspect registered export specs
     (`MAPL_Get(STATE, ExportSpec=...)`), so a scoped
     `MAPL_GetObjectFromGC` was added there with a `TODO (step 4)`
     comment. Step 4 removes it along with the one in `SetServices`.

8. **Edge remaps** - `PLE`/`ZLE` are `VLOC=E`; MAPL3 pointers are
   1-based but `SIM_DRIVER` indexes `PLE(:,:,0:LM)` and declares
   `PLE2D(IM*JM,0:LM)` locals. Apply the contiguous 0-based rank-remap
   trick from Solar step 11 / IRRAD notes, or shift the indexing.
   - DONE 2026-10-22. `PLE` is remapped right after
     `#include "Satsim_GetPointer___.h"` in `SIM_DRIVER`:
     `p3d => PLE; PLE(1:IM, 1:JM, 0:LM) => p3d`, with a local
     `real, pointer, contiguous, dimension(:, :, :) :: p3d` (a self-remap
     does not compile, IRRAD notes). No `associated()` guard: `PLE` is a
     mandatory import. Its 0-based uses are `PLO`/`DELP`
     (`PLE(:, :, 0:LM - 1)`) and `PLECOSP` (`PLE(:, :, LM:0:-1)`).
   - `ZLE` is remapped the same way, for consistency: it is only used as
     `reshape(ZLE, ...)` into the 0-based local `ZLE2D`, which does not
     depend on its bounds, but the remap keeps both edge imports 0-based
     should anyone index `ZLE` directly. These are SatSim's only edge
     fields (all 222 exports are `C` or `N`).
   - The ACG does NOT declare its pointers `contiguous` (checked
     `MAPL_GridCompSpecs_ACG.py` and the generated
     `Irrad_DeclarePointer___.h`); Solar step 11's correction is right,
     IRRAD's "fixed at the ACG root" note is out of date. The scratch
     pointer is what makes the remap legal.
   - Not compiled yet (step 11).

9. **Entry-point signatures** DONE 2026-10-01.
   - `SetServices(GC, RC)` now takes a plain `type(ESMF_GridComp) :: GC`
     and `integer, intent(out) :: RC`, matching MAPL3's `I_SetServices`
     (`infrastructure/esmf/ESMF_Interfaces.F90`). The MAPL2 version had
     `intent(inout) :: GC` and an `optional :: RC`, neither of which is
     allowed by the abstract interface.
   - Added an `Initialize(GC, IMPORT, EXPORT, CLOCK, RC)` with the
     `I_Run` signature (all dummies plain, `rc` `intent(out)`); `Run`
     was changed to match. Both are registered with
     `MAPL_GridCompSetEntryPoint(GC, ESMF_METHOD_INITIALIZE/RUN, ...)`,
     the Run one with `phase_name="run"` as for IRRAD/SOLAR. Finalize
     uses the MAPL generic default.
   - Dropped the trailing `MAPL_GenericSetServices(GC, _RC)`; MAPL3
     wires the generic methods itself.
   - `Initialize` creates the run alarm (step 7's alarm half, done here
     since the new method had to exist anyway): `SATSIM_DT` seconds,
     defaulting to the heartbeat, with `_ASSERT`s that it is positive
     and a multiple of the clock step (a zero interval would hang the
     `ring_time` loop). The ring-time arithmetic is copied verbatim from
     IRRAD/SOLAR, i.e. MAPL2's `handle_clock_and_main_alarm`, honouring
     `RUN_AT_INTERVAL_START` / `REFERENCE_DATE` / `REFERENCE_TIME`. The
     alarm is named `satsim_alarm`; step 7 will have `Run` fetch it with
     `ESMF_ClockGetAlarm` in place of the old generic `RUNALARM`.
   - `_UNUSED_DUMMY(IMPORT)` / `_UNUSED_DUMMY(EXPORT)` in `Initialize`,
     since the states are untouched there.
   - The container's `use ... only: satsimSetServices => SetServices`
     alias still suffices (step 10).

10. **Build wiring** - DONE 2026-10-22.
    - `GEOSsatsim_GridComp/CMakeLists.txt`: no change. `DEPENDENCIES`
      already lists `GEOS_Shared MAPL GEOS_RadiationShared`. The plan's
      "`TYPE SHARED` to match siblings" was wrong: IRRAD and Solar don't
      set it (only the parent does).
    - DONE - parent `CMakeLists.txt`: `GEOSsatsim_GridComp` uncommented
      in `alldirs`.
    - DONE - `../GEOS_RadiationGridComp.F90`:
      - `use GEOS_SatsimGridCompMod, only: satsimSetServices =>
        SetServices` restored.
      - `SetServices`: as in MAPL2, `USE_SATSIM` is the sum of the six
        `USE_SATSIM*` flags (each default 0, read by `get_use_satsim`),
        and `USE_SATSIM > 0` gates `MAPL_GridCompAddChild(gc, "SATSIM",
        satsimSetServices, "satsim.yaml", _RC)`. The flags only need to
        be set once, in the parent's config: MAPL3 merges the parent's
        non-`mapl:` keys into the child's config, with the child's value
        winning if both define a key (`add_child_by_spec.F90:38`,
        `merge_hconfig` in `HConfigUtilities.F90:24-59`), so SATSIM sees
        the same `USE_SATSIM*` values.
      - One thing to watch: when parent and child both define a key, the
        child's value wins (`HConfigUtilities.F90:51-52`). A `satsim.yaml`
        used under RADIATION should therefore not set the `USE_SATSIM*`
        flags. If it did, SATSIM would use its own values, which could
        differ from the ones RADIATION used to decide whether to create
        SATSIM.
      - `MAPL_GridCompAddConnection(gc, src_comp="SOLAR",
        src_names="MCOSZ", dst_comp="SATSIM")`, MAPL2's only explicit
        SATSIM connection (Solar exports `MCOSZ` with the same
        standard name). The other 18 imports propagate up unconnected,
        as IRRAD's and Solar's do.
      - `Run`: SATSIM runs after SOLAR (MAPL2 order SOL, IRR, STM; it
        needs this step's `MCOSZ`), behind the same `USE_SATSIM > 0` check.
      - Timers fixed in the parent's `Run` while there: the port had
        turned MAPL2's "pause `TOTAL` around the children"
        (`MAPL_TimerOff(TOTAL)`, `MAPL_GenericRunChildren`,
        `MAPL_TimerOn(TOTAL)`) into per-child timers in inverted order
        (`TimerStop("IRRAD")`, `RunChild`, `TimerStart("IRRAD")`), and
        dropped the `TimerStart("TOTAL")` at the top, leaving only the
        `TimerStop("TOTAL")` at the end. Restored MAPL2's pattern: start
        `TOTAL` on entry, stop it around the three `RunChild` calls, stop
        it at the end; no per-child timers.
    - No re-exports: MAPL2 never promoted SatSim exports into RADIATION
      (no `CHILD_ID=STM`); they live under `SATSIM`, which is how the
      `HISTORY*.rc.tmpl` templates request them, and MAPL3 propagates
      child exports under the child's name. The plan's "re-export the
      224+ exports" was wrong.
    - `satsim.yaml` (`SATSIM_DT`, `USE_SATSIM_*`, `SATSIM_NCOLUMNS`) only
      exists where a run needs it; first one in step 12. `SATORB` and
      `Masked_Exports` are not needed while masked exports are disabled
      (step 4).
    - Decided 2026-10-22: SatSim's import metadata is kept as in MAPL2,
      NOT aligned with IRRAD/Solar. It differs for 15 shared imports
      (13 standard names, 4 units), e.g. `PLE` `Edge_pressures` vs
      `air_pressure`, `TS` `skin_temperature` vs
      `surface_skin_temperature`, `QV`/`QL`/`QI` units `1` vs `kg kg-1`,
      `RR` `rain_cloud_particle_effective_radius` vs
      `effective_radius_of_rain_particles`. These imports are merged
      under RADIATION and connected to one provider, so expect a
      standard_name warning per mismatch in the default permissive mode
      (the provider's value wins, and is what checkpoints record) and
      errors in strict mode. Matters for step 12: the `satsim-sa`
      import checkpoint will carry the fake provider's names.
    - 2026-10-26, `regression/radiation` runs SatSim (ISCCP only, as
      `satsim-sa`). Its `DataForRadiation` declares the shared imports
      with IRRAD/Solar metadata, plus SatSim-only `ZLE`, `FRLAND` and
      `FROCEAN` with SatSim's, each connected to RADIATION. `MCOSZ`
      comes from Solar. As expected above, SatSim's shared imports log
      standard-name warnings (`StandardNameAspect.F90:150-163`,
      permissive is the default, `FieldDictionaryConfig.F90:92`), and
      `QV`/`QL`/`QI` get a `kg kg-1` -> `1` units conversion
      (`UnitsAspect.F90:89-91`) that leaves the values unchanged: all
      18 provider-fed imports in SatSim's checkpoint equal the input.
      `satsim-sa` now nests SATSIM under a stand-in RADIATION
      (configurable gridcomp) that declares the `USE_SATSIM*` flags, and
      `satsim.yaml` no longer sets them (see the "thing to watch" above).

11. **Compile with ifx** and fix fallout. Only the GridComp should
    break; the fixed-form COSP `.f` sources are untouched.
    - DONE 2026-10-22: the build with SatSim wired in (step 10) succeeded
      with no code changes needed (built by the user). Not run yet
      (step 12).

12. **Standalone regression `regression/satsim-sa`** - clone
    `irrad-sa` (`cap_*.yaml`, `mapl.yaml`, `root.yaml`, `data-*.yaml`
    provider, `cmpchk.py`). The fake provider must supply all 19 field
    imports (no `SATORB` bundle while masked exports are disabled, step 4), using dictionary
    `standard_name`s (Solar notes RR lesson). Generate a MAPL2 baseline
    (SatSim has never been run-tested in this repo, same as IRRAD was)
    and add the case to `regression/cases.txt` for ctest.
    - The MAPL2 baseline must match the standalone's configuration or it
      will never reproduce bit-for-bit (Solar step 18, SU step 0): record
      every `USE_SATSIM*`, `SATSIM_NCOLUMNS`,
      `SATSIM_POINTS_PER_ITERATION` and `SATSIM_DT` used, and request
      the compared exports in `HISTORY.rc` (MAPL2 allocates exports only
      on request).
    - Exports a run never writes are NaN in MAPL3 where MAPL2 had zeros
      or stale values (SU step 9); expect that for simulators that are
      switched off.
    - SatSim's exports use `MAPL2_UNDEF = 1e15`, as in MAPL2 (step 4);
      `cmpchk.py` should still treat undef cells as equal.
    - Hand runs need `LD_LIBRARY_PATH` set to the install `lib`, and
      `cap_restart.yaml` reset before each run (Solar step 18).
    - DONE 2026-10-22 - case files in `../regression/satsim-sa/`, after
      `solar-sa` (not run yet, no baseline yet, not in `cases.txt`):
      - `mapl.yaml`, `cap_driver.yaml` (PT20M), `cap_gridcomp.yaml`,
        `cap_restart.yaml` (2000-04-14T22:00) copied from `solar-sa`; the
        baseline must be captured at the matching time.
      - `root.yaml`: C12L181 geometry, with the current vertical-grid
        schema (`grid_type: model`, `fields: {pressure: PLE}`, as in
        `GEOSphysics_GridComp/regression/radiation+gocart`. MAPL no
        longer reads the old `class: model`/`field_edge` form; on
        2026-10-22 `irrad-sa`, `solar-sa` and FV's `dyn-sa`/`adv-dyn`
        were converted too, and all regression tests pass); children
        `DataForSatsim` (configurable gridcomp) and `SATSIM`
        (`libGEOSradiation_GridComp.so`, `satsim_setservices_`); one
        connection per import from `DataForSatsim`'s internal state. In
        the full model `MCOSZ` comes from SOLAR; here the provider has
        it too.
      - `data-for-satsim.yaml`: the 19 imports as restarted internals,
        with SatSim's own units and standard names (step 10 decision),
        so they match the MAPL2 baseline's import metadata. Its restart
        (`checkpoints/2000-04-14T22:00:00/DataForSatsim_internal.nc`)
        comes from the MAPL2 SATSIM import state before the run.
      - `satsim.yaml`: `SATSIM_DT: 1200`, `USE_SATSIM: 0`,
        `USE_SATSIM_ISCCP: 1` (ISCCP only; the MAPL2 capture run cannot
        use `USE_SATSIM: 1`), other per-simulator flags 0,
        `SATSIM_NCOLUMNS` and `SATSIM_POINTS_PER_ITERATION` left at their
        defaults (30 at C12, -999); `activate_all_exports: true`, import
        and export checkpoints (SatSim has no internal state). The MAPL2
        capture must use the same values.
      - MAPL2 `HISTORY.rc` may request only the exports MAPL2 fills with
        ISCCP alone: the 69 ISCCP exports and `SGFCLD` (70 of 222). Each
        simulator's local arrays are copied out of COSP only inside that
        simulator's `USE_SATSIM + USE_SATSIM_<SIM> > 0` block
        (MAPL2 `GEOS_SatsimGridComp.F90` L3738-3790), but the
        `if (associated(X))` export assignments after them are ungated.
        So requesting e.g. `LIDARPTOT` reads the never-written
        `lidar_beta_tot`, which a Debug build (`-init=snan,arrays
        -fpe0`) traps as "floating invalid" at L4330. `SGFCLD` is safe:
        `sgx` is built unconditionally (L3647). Export collections:
        `geosgcm_satsim` (56 fields) and `geosgcm_satsim_ungrd` (the 14
        size-7 `CLISCCP1-7`, `ISCCP1-7`). A MAPL2 HISTORY collection with
        no vertical regridding requires every rank-3 field to have the
        same 3rd-dim size, ungridded dims included (MAPL
        `base/MAPL_VerticalMethods.F90:666-669,713-724,752`, "have mixed
        level/edge"), so ungridded fields go in their own collection.
      - OPEN - the MAPL3 port has the same ungated export assignments,
        and `activate_all_exports: true` allocates every export, so with
        ISCCP alone 107 exports are computed from never-initialized work
        arrays: `LIDARPTOT`/`LIDARTAUTOT` (`lidar_beta_tot`,
        `lidar_tau_tot`), the 34 MISR and the 71 MODIS exports. The
        "zero out exports" block (`GEOS_SatsimGridComp.F90:717-734`)
        covers only the ISCCP, `frac_out`, radar and other lidar arrays.
        A Debug build (`-init=snan,arrays -fpe0`,
        `@cmake/compiler/flags/IntelLLVM_Fortran.cmake:187`) traps at
        `LIDARPTOT` (L1399) as MAPL2 did; a Release build writes
        garbage, which breaks `run_case`'s whole-file
        `nccmp -dmfgsB --nans-are-equal`
        (`@cmake/esma_support/esma_regression_run_helpers.cmake:92-121`).
        `cmpchk.py` itself is unaffected: it compares only variables
        present in both files.
      - `cmpchk.py`: Solar's, for the `import` and `export` states of
        `SATSIM`, against `~/input/satsim/C12L181/after/new-style`.
    - DONE 2026-10-22 - added the external `Satsim_SetServices(gc, rc)`
      wrapper after `end module`, as Solar step 13, so the standalone can
      load SatSim by symbol (`satsim_setservices_`). Not compiled yet.
    - No ExtData needed: `gettau` only uses constants, and the COSP file
      reads (`icarus.f` debug output, `geos5_read.f90`) are not on the
      driver's path.

13. **Style cleanup afterwards** (low priority, like Solar step 16).
    DONE 2026-10-01 (done out of order, before steps 3-12).
    - Keywords/alignment: nothing to do, the step-1 `codee format` pass
      had already lowercased every Fortran keyword and collapsed all
      column-aligned `::` declarations. Only 7 uppercase *intrinsics*
      remained; lowercased them (`AdjustL`, `MIN`/`MAX`/`INT` in the
      `default_Ncolumns` associate, 4 `MAX` in the `CWC = max(Q/FCLD,
      1e-12)` block, `EXP` in `EMISS = 1 - exp(-taucir)`).
    - Removed the commented-out `frac_outinv` declaration/allocate/
      deallocate (3 lines), the commented `construct_cosp_sghydro` /
      `FREE_COSP_SGHYDRO` pair plus the now-unused
      `type(cosp_sghydro) :: sghydro` declaration, the 5 commented
      debug `write(*,*)` lines inside the `DEBUG_GC` block, and the
      13-line commented-out pre-`SATSIMTEMP` `RADARZETOT` block.
    - `CMakeLists.txt`: dropped the commented-out
      `#  quickbeam/load_mie_table.f90` and `#  congvec.f` `srcs`
      entries. **`congvec.f` itself must stay on disk** - `scops.f`
      lines 177 and 246 still `include 'congvec.f'`, so it is a text
      include, not a compiled unit. `load_mie_table.f90` is genuinely
      unreferenced (only a comment in `cosp_types.F90:1134`).
    - Fixed the `nullify`-inside-loop bug in the masking pass (the
      `MAPL_DimsHorzOnly .and. associated(ungridded_dims)` branch
      nullified `ptr3d_new`/`ptr3d` inside the
      `do j = 1, ungridded_dims(1)` loop): moved the two `nullify`
      calls out of the loop, matching the other branches. Done now
      rather than "once a test covers it" - the regression test is
      still step 12, so this one is reviewed-by-inspection.
    - Left alone deliberately: `DEBUG_GC` and the live `write(*,*)`
      debug prints (step 7 removes those) and the `! deallocate stuff`
      section comment.
    - Not recompiled at the time (step 11 was still open then).
    - 2026-10-22 follow-up (after the audit). Builds, and all
      regression tests pass. SatSim itself has not been run
      (`satsim-sa` is not a regression case yet, step 12):
      - `#include "MAPL_Generic.h"` -> `"MAPL.h"`, as IRRAD/RADIATION.
        `MAPL_Generic.h` is only `#include "MAPL.h"`
        (`@MAPL/include/MAPL_Generic.h:1`), so no behavior change.
      - Removed unused locals, each confirmed by a per-routine use
        count over the source plus the expanded ACG headers, and a
        separate line-by-line search: `DIMS` (`SetServices`); `name`,
        `frac_ls`, `l`, `scops_debug`, `idebug`, `idebugcol`,
        `TAUDIAG`, `TAUSWICE`, `TAUSWLIQ`, `EMISS2D` (`SIM_DRIVER`).
        The `SatSim_State` components were left alone (private state
        kept on purpose, step 5).
      - `_RC` on the one call that takes `rc` and lacked it:
        `alarm_ringing = ESMF_AlarmIsRinging(satsim_alarm, _RC)` before
        the early return in `Run`, as IRRAD does
        (`GEOS_IrradGridComp.F90:680`). `_RC` expands to
        `rc=status);_VERIFY(status` (`MAPL_ErrLog.h:114`), so it
        cannot sit inside an `if (...)` condition; hence the local.
        `MAPL_UngriddedDim` (pure) takes no `rc`.
      - `Initialize` stays registered without `phase_name`: MAPL runs
        the user Initialize by the default name `GENERIC::INIT_USER`
        (`initialize_user.F90:21`), and `run_custom` skips a phase it
        can't find without error (`_RETURN_UNLESS(found)`), so a custom
        name would silently disable it. All other ports do the same.

14. **SatSim writes into its radius imports** - DONE 2026-10-26 (found
    in the 2026-10-22 audit; the original notes follow, then the
    resolution).
    - `SIM_DRIVER` (~L543-557) assigns to the import pointers
      themselves: `where (RDFI == MAPL2_UNDEF) RDFI = 36.e-6`, and the
      same for `RDFL`, `RDFR`, `RDFS`, `RDFG` (imports `RI`, `RL`,
      `RR`, `RS`, `RG`). Unchanged from MAPL2.
    - Until the `MAPL2_UNDEF` fix (step 4) these tests never matched,
      so the writes were inert in MAPL3. They now fire wherever the
      inputs hold `1e15`.
    - Under RADIATION, IRRAD and Solar import the same fields. SatSim
      runs after both, so the current step is not affected.
    - Unverified, to check:
      - whether a MAPL3 connected import shares memory with the
        provider's export (MOIST's `RI`... in the full model), in
        which case the replaced values reach every other consumer,
        HISTORY, and the next step's IRRAD/Solar;
      - whether MAPL2 shared that memory too (if so, MAPL2 had the
        same side effect and the baseline includes it);
      - whether anything downstream depends on `1e15` marking
        cloud-free radii.
    - Likely fix: copy the radii into locals before replacing the
      undef values, as Solar does (`RR3`/`REFF` are locals).
    - Resolved 2026-10-26:
      - Confirmed by the `radiation` regression case once SATSIM was
        enabled in it (ISCCP only): the end-of-step `IRRAD_import.nc`
        checkpoint had SatSim's defaults in place of `1e15`, at all
        121,436 `RI` cells (36e-6) and 132,799 `RL` cells (14e-6).
        So IRRAD's, Solar's and SatSim's imports share one field under
        RADIATION (no units conversion on the radii, all `m`), and
        SatSim's writes reach every other reader after it runs.
        IRRAD/Solar exports were unaffected, since they run first.
        `RR`/`RS`/`RG` hold no `1e15` in that input, so did not change.
      - Only these five imports are written. Nothing else in SatSim
        assigns to an import, and the only imports passed to
        subroutines, `MCOSZ` and `FCLD` (`getvistau`/`getirtau`), are
        `intent(IN)` (`../GEOS_RadiationShared/gettau.F90:45,47,187`).
      - MAPL2 writes into the import pointers too (`origin/develop`
        `GEOS_SatsimGridComp.F90:3080-3084`, `:3394-3410`). Whether
        MAPL2 also shares that memory with IRRAD/Solar is still
        unverified. The MAPL2 capture cannot tell: each component's
        "after" state was saved right after it ran, before SatSim, so
        IRRAD's and Solar's `RI`/`RL` hold `1e15` while SatSim's import
        holds the defaults. Reported for MAPL2 as
        [#111](https://github.com/GEOS-ESM/GEOSradiation_GridComp/issues/111).
      - Fix: `SIM_DRIVER` copies the five imports into locals
        `RDFI_LOCAL` ... `RDFG_LOCAL` (declared ~L314-316) and applies
        the defaults to those (~L538-558). `REFF` (~L561-564) and the
        `RDF*COSP` arrays (~L699-703) read the locals. The imports are
        now only read.
      - Effect: exports unchanged (bit-for-bit in `satsim-sa`). The
        SatSim import checkpoint now keeps `1e15` in the radii, so it
        no longer matches the MAPL2 capture's SatSim import in
        `RI`/`RL`, and `cmpchk.py` reports ~1e15 differences there
        (it treats only NaN, `+-huge` and `_FillValue` as undef).
      - Regression data (`~/regression-data`, 2026-10-26; `/ford1` not
        updated yet): `satsim-sa/checkpoints/last/SATSIM_import.nc`
        replaced (the old one differed only in `RI`/`RL`).
        `radiation/checkpoints/last` gained `SATSIM_import.nc` and
        `SATSIM_export.nc`, from a hand run that matched all eight
        existing files. The `radiation` SatSim export and import are
        identical to `satsim-sa`'s, including `MCOSZ`, which comes from
        Solar there.

Effort concentrates in steps 2 (mechanical but ~2300 lines), 4 (needs a
MAPL3 design decision for dynamic export cloning) and 12 (test data and
baseline).
