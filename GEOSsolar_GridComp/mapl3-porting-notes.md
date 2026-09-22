# GEOS_SolarGridComp MAPL3 porting plan

Status: in progress. `GEOSsolar_GridComp` is now uncommented in the
parent container's `alldirs` (see `../CMakeLists.txt`) and wired as a
child of `GEOS_RadiationGridComp.F90`. Steps 1-11 and 13-15 below are
done (see branch `feature/pchakrab/port-solar-to-mapl3`); step 12 has
been investigated but not yet implemented (see its notes below). The
port has **not yet been build-verified end-to-end**.

This plan was derived by comparing against the already-completed IRRAD
port (`../GEOSirrad_GridComp/`, see its own `mapl3-porting-notes.md` for
the full list of MAPL3 API discoveries/gotchas found while doing that
port - referenced throughout below instead of repeated).

## Reference files
- Template component: `../GEOSirrad_GridComp/GEOS_IrradGridComp.F90`
- Template state spec: `../GEOSirrad_GridComp/Irrad_StateSpecs.rc`
- Template build file: `../GEOSirrad_GridComp/CMakeLists.txt`
- Full IRRAD lessons-learned log: `../GEOSirrad_GridComp/mapl3-porting-notes.md`
- Parent container (already updated for IRRAD, needs Solar wiring added):
  `../GEOS_RadiationGridComp.F90`
- Parent's own porting notes (documents exactly what's deferred until
  Solar exists): `../mapl3-porting-notes.md`

## ACG script notes (from step-3 investigation)
- Two different ACG scripts exist in this environment - don't confuse
  them. The only one that matters for final validation of
  `Solar_StateSpecs.rc` is the real, in-repo
  `../../../../../../../Shared/@MAPL/apps/MAPL_GridCompSpecs_ACG.py`
  (i.e. `src/Shared/@MAPL/apps/MAPL_GridCompSpecs_ACG.py` from the repo
  root). It requires Python 3.10+ (`match` statement) - use
  `python3.12`, not the system default `python3` (3.6.8, too old).
  CLI: `python3.12 MAPL_GridCompSpecs_ACG.py Solar_StateSpecs.rc -i
  imp.h -x exp.h -p int.h` (no `-n`/`--debug` flags). A separate
  MAPL2-era reference copy at `~/mapl2/solar/MAPL_GridCompSpecs_ACG.py`
  is useful only for understanding legacy semantics, not for final
  validation - its schema is meaningfully different (see below).
- MAPL3's ACG generates one unified
  `call MAPL_GridCompAddSpec(gridcomp=gc, short_name=, units=, dims=,
  vertical_stagger=MAPL_VERTICAL_STAGGER_{CENTER,EDGE,NONE},
  standard_name=, state_intent=ESMF_STATEINTENT_{IMPORT,EXPORT,
  INTERNAL}, _RC)` call per spec, not separate `MAPL_AddImportSpec`/
  `AddExportSpec`/`AddInternalSpec` calls like MAPL2. The `.rc`'s
  `LONG NAME`/`LONG_NAME` column feeds `standard_name=` (not a separate
  `long_name=` arg - a schema rename vs MAPL2).
- Columns confirmed UNSUPPORTED by MAPL3's ACG (silently dropped, no
  error/warning - `main()` only reports missing *mandatory* keys, never
  unrecognized ones): `FRIENDLYTO`, `DEFAULT`, `AVERAGING_INTERVAL`,
  `REFRESH_INTERVAL`, `ACCMLT_INTERVAL`, `COUPLE_INTERVAL`.
- `FRIENDLYTO`/`DEFAULT` and `AVERAGING_INTERVAL`/`REFRESH_INTERVAL` get
  *different* treatment despite both being unsupported columns:
  - **`FRIENDLYTO` (CORRECTED)**: `FRIENDLYTO=trim(comp_name)` on a MAPL2
    INTERNAL spec is NOT inert - it's the mechanism (per Solar's own
    comment at `GEOS_SolarGridComp.F90:742-747`) that auto-exports an
    INTERNAL under the *same* short_name with zero explicit EXPORT
    declaration and zero manual copy code. MAPL3's real functional
    replacement is the INTERNAL-category **`ADD2EXPORT`** column ->
    `add_to_export=.true.` arg on `MAPL_GridCompAddSpec`
    (`gridcomp_add_spec` in `MAPL_Generic.F90:622-758`). When
    `add_to_export=.true.` and `export_name` is omitted,
    `gridcomp_reexport`/`ComponentSpec%reexport`
    (`specs/ComponentSpec.F90:122-158`) defaults the new export name to
    `src_name` - i.e. auto-export under the same name, exactly matching
    MAPL2 FRIENDLYTO semantics. Confirmed real usage:
    `Irrad_StateSpecs.rc`'s INTERNAL category has an `ADD2EXPORT` column;
    its `DFDTS` row has `ADD2EXPORT=.TRUE.` with a comment "DFDTS is
    declared in the INTERNAL category above with ADD2EXPORT so it is
    automatically also exported under the same name" - and `DFDTS` does
    NOT appear as a separate EXPORT row. **Rule for Solar: every INTERNAL
    with `FRIENDLYTO=trim(comp_name)` that has no separate EXPORT row of
    the same name should get `ADD2EXPORT=.TRUE.` in `Solar_StateSpecs.rc`,
    not be treated as inert/dropped.** (Solar has ~149 FRIENDLYTO
    occurrences - this affects a large fraction of the INTERNAL category.)
    `GEOS_Ocean_StateSpecs.rc`'s `TS_FOUND` row (which just keeps
    `FRIENDLYTO=trim(COMP_NAME)` as a literal, ACG-ignored column with no
    `ADD2EXPORT`) looks like an incomplete/lazy port and should NOT be
    followed as precedent.
  - `DEFAULT`: the underlying `MAPL_FriendlyVariable`-adjacent default
    mechanism appears fully dropped in MAPL3 ACG (no real replacement
    column found so far); `GEOS_Ocean_StateSpecs.rc` keeps it as inert
    documentation and the real ACG accepts this with zero errors. Treat
    as inert/droppable unless a functional replacement turns up.
  - `AVERAGING_INTERVAL`/`REFRESH_INTERVAL` (MAPL2 IRRAD/Solar's
    `ACCUMINT`/`MY_STEP`): confirmed via the MAPL2-era ACG-fy of IRRAD
    on `origin/refactor/pchakrab/acgfication` (in the
    `@GEOSradiation_GridComp` nested repo; pre-ACG-fy commit `5241d6a`
    has the manual `MAPL_AddImportSpec(..., AVERAGING_INTERVAL=
    ACCUMINT, REFRESH_INTERVAL=MY_STEP, _RC)` calls, ACG-fy commit
    `79993dc` "ACG-fy SetServices state specs" generated a MAPL2-schema
    `.rc` that DID preserve these as real columns) that this was a
    valid MAPL2 ACG feature. But the actual MAPL3 port of IRRAD in
    this repo has ZERO trace of AVERAGING_INTERVAL/REFRESH_INTERVAL/
    ACCUMINT/MY_STEP anywhere (not even as an inert `.rc` column) -
    full removal, apparently superseded by ExtData2G-driven time-
    accumulation. **Established precedent for Solar: drop
    AVERAGING_INTERVAL/REFRESH_INTERVAL/ACCMLT_INTERVAL/
    COUPLE_INTERVAL entirely** (no `.rc` column, no `ACCUMINT`/
    `MY_STEP`-equivalent variables/logic in `SetServices`) - follow
    IRRAD's full-removal treatment, not Ocean's
    inert-documentation-column treatment.

## Steps

1. **Reformat for consistent indentation and RC macro style, before any
   functional changes.** Do this first so later diffs (state-spec
   extraction, signature fixes, etc.) are clean and don't mix
   whitespace/style noise with real changes. Matches the order IRRAD's
   own port actually happened in (see its notes file's "whitespace
   cleanup" session, done early alongside the `SetServices` signature
   fix, well before the bulk of the functional MAPL3 conversion).
   - `__STAT__` -> `_STAT`, `__RC__`/`VERIFY_`/`RETURN_`/`ASSERT_` ->
     `_RC`/`_VERIFY`/`_RETURN`/`_ASSERT` throughout (both are
     functionally identical under this build's macro definitions since
     ANSI_CPP is not defined - pure naming-convention change, verify
     old/new macro expansions match before bulk-replacing).
   - Fixed continuation-line indentation: this file's (and IRRAD's)
     established convention is a **fixed 5-space offset** from the
     starting statement's own indent for every `&`-continuation line -
     not alignment to the opening parenthesis. Watch for continuation
     lines with a trailing inline comment after the `&`
     (`.false. , &! comment`) - a naive "ends with `&`" check
     under-detects those; strip trailing `!...` comments first
     (tracking quote state) before testing `endswith('&')`.
   - Lowercase stray ALL-CAPS keywords (`IF`/`THEN`/`ENDIF`/`WHERE`/
     `ALLOCATE`/`DEALLOCATE`/`REAL`/`INTEGER`/`POINTER`/`DIMENSION`/
     etc.) and `intent(IN )`/`intent(OUT)` -> `intent(in)`/`intent(out)`,
     matching the lowercase-keyword convention used file-wide in IRRAD.
     Leave OpenMP sentinel comments (`!$OMP ...`) and non-keyword
     identifiers (e.g. resource-label names) alone.
   - Collapse column-alignment padding in `use ... only:` statements and
     declaration blocks (interior multi-space runs used to vertically
     align `intent(...)`/`only:`/closing parens across neighboring
     lines) down to single-space.
   - Remove pure dash-separator comment lines (`!----------` style)
     that are purely decorative.
   - Verify every one of the above is whitespace/token-rename-only via
     a strip-and-diff (or tokenize-and-compare-identifier-multiset for
     declaration-block reordering/merging) pass - no logic change
     should occur in this step.
   - Do NOT attempt this as a single whole-file automated fprettify run
     - IRRAD's notes document that fprettify re-flows already-correct
     continuation styles it doesn't understand (paren-alignment vs. the
     fixed-5-space convention used here). If using fprettify at all,
     scope each invocation to one subroutine at a time and diff
     before/after to confirm the change is pure whitespace.

2. **DONE - `SORADCORE`/`UPDATE_EXPORT` contained-procedure scope.**
   Like IRRAD's `LW_Driver`, these two are large procedures nested in
   `Run`'s `contains`, referencing dozens of host-associated pointers -
   left nested (per IRRAD's lesson, don't hoist these with a giant
   explicit-argument list). The RRTMGP-only helper cluster
   (`compute_gas_optics`, `compute_aer_optics`,
   `compute_cloud_optics_mcica`, `compute_sprlyr_diags_predelta`,
   `compute_delta_scale`, `compute_sprlyr_diags_postdelta`,
   `compute_rte_sw`, `PROCESS_RRTMGP_BLOCK`, `shrtwave`) was verified
   self-contained and has been hoisted to module scope, placed
   immediately after `end subroutine Run`. `shrtwave` needed three new
   explicit dummy arguments (`HK_UV_TEMP`, `HK_IR_TEMP`, `MAPL`) that
   had been host-associated from `Run`'s locals; its sole call site
   (inside `SORADCORE`) was updated to pass them explicitly. Verified
   via before/after structural marker counts (subroutine/end-subroutine,
   `#ifdef`/`#endif`, `#define`/`#undef TEST_` pairs) and grep for
   residual host-association - Solar isn't wired into the build yet so
   this could not be compile-verified.

2a. **DONE - codee format pass over the step-2 changes.** Ran `codee
    format` on the whole file to bring the moved/modified code in line
    with the file's style (see
    `.vscode/instructions/fortran-formatting.instructions.md` for the
    mask/format/unmask workaround needed for a handful of constructs
    codee's parser can't handle). Verified via the same before/after
    structural marker counts as step 2. Landed together with step 2 as
    a single commit on `feature/pchakrab/port-solar-to-mapl3`.

3. **DONE - Extract state specs into `Solar_StateSpecs.rc`.** Pulled every
   `MAPL_AddImportSpec`/`AddExportSpec`/`AddInternalSpec` call out of
   `GEOS_SolarGridComp.F90`'s `SetServices` into
   [Solar_StateSpecs.rc](Solar_StateSpecs.rc) (`IMPORT`/`EXPORT`/
   `INTERNAL` categories). The `AERO` nested state (an `ESMF_State`, not
   a plain Field - ACG's `ITEMTYPE` column only supports `F`/`V`) and the
   `OSRBbbRG`/`ISRBbbRG`/`TBRBbbRG` per-band dynamic specs (loop over
   `ibnd`) are genuinely dynamic/non-tabular - left as manual
   `MAPL_GridCompAddSpec` calls in `SetServices`, same treatment as
   IRRAD's RATS diagnostics loop.

4. **DONE - Add `mapl_acg()` to `CMakeLists.txt`** with `IMPORT_SPECS
   EXPORT_SPECS INTERNAL_SPECS GET_POINTERS DECLARE_POINTERS`, and added
   `TYPE SHARED` to `esma_add_library()`.

5. **DONE - Replace static `MAPL_Add*Spec` blocks** in `SetServices` with
   `#include "Solar_Import___.h"` / `_Export___.h` / `_Internal___.h`.

6. **DONE - Convert the RRTMGP internal state.** Replaced the manual
   `ty_RRTMGP_wrap`/`ESMF_UserCompSetInternalState`/
   `ESMF_UserCompGetInternalState` pattern with `_SET_NAMED_PRIVATE_STATE`/
   `_GET_NAMED_PRIVATE_STATE(gc, ty_RRTMGP_state, PRIVATE_STATE, ...)`
   macros at all 4 call sites (`SetServices`, and 3 sites inside `Run`'s
   contained procedures `SORADCORE`/`UPDATE_EXPORT`, which are
   host-associated with `gc` so no signature changes were needed). Added
   a module-scope `character(*), parameter :: PRIVATE_STATE =
   "RRTMGP_state"` and removed the now-dead `ty_RRTMGP_wrap` type,
   matching IRRAD's pattern exactly.

7. **DONE - `MAPL_GetResource(MAPL,...)` -> `MAPL_GridCompGetResource(gc,
   "LABEL", var, default=..., _RC)`** (drop trailing colon on labels)
   throughout `SetServices` and `Run`. Converted all 40 call sites
   (`SetServices`, `Run`, and inside `SORADCORE`/`UPDATE_EXPORT`/
   `PROCESS_RRTMGP_BLOCK`, all host-associated with `gc`).

8. **DONE - `ESMF_Config` -> `ESMF_HConfig`.** Solar never read config
   directly via `ESMF_Config` (only through `MAPL_GetResource`/now
   `MAPL_GridCompGetResource`) - removed the now-unused `type(ESMF_Config)
   :: cf` declaration and the `cf=cf` argument to `MAPL_Get` in `Run`;
   no `ESMF_HConfig` usage needed.

9. **DONE - Fix entry-point signatures**: `SetServices(gc, rc)`/
   `Run(gc, import, export, clock, rc)` with **no `intent`/`optional`**
   on `gc`, and `rc` non-optional `intent(out)`. This is required
   before `MAPL_GridCompAddChild` can add Solar as a child from
   `GEOS_RadiationGridComp.F90`. Dropped `intent(inout)` from `gc` and
   `optional`/`intent(inout)` from the other dummy args in both
   `SetServices` and `Run` to match `Irrad`'s signatures exactly;
   confirmed no `present(rc)` checks existed anywhere in the file, so
   making `rc` non-optional is safe.

10. **DONE - Convert `MAPL_MetaComp`/`MAPL_Get`/`MAPL_TimerOn` usage in
   `Run`.** `MAPL_GenericSetServices` (the MAPL2 mechanism that let a
   component define only `Run` and get a default `Initialize`/
   `Finalize`) is completely gone from MAPL3 - confirmed zero matches
   anywhere in `src/Shared/@MAPL`, including the Deprecated shim. Since
   `SetServices` already registered `ESMF_METHOD_INITIALIZE, Initialize`
   (a leftover from before this port that referenced a since-removed
   subroutine - a latent bug), a **new `Initialize(gc, import, export,
   clock, rc)` subroutine had to be added** (mirroring IRRAD's), which:
   - Creates a `"solar_run_alarm"` (`ESMF_AlarmCreate`, `sticky=.true.`)
     with ring interval from `<NAME>_DT` (default: heartbeat via
     `MAPL_ClockGet(clock, dt=..., _RC)`), matching IRRAD's
     `"irrad_lw_alarm"` pattern. `Run` retrieves it via
     `ESMF_ClockGetAlarm(clock, alarmname="solar_run_alarm", alarm=alarm, _RC)`
     in place of MAPL2's `MAPL_Get(MAPL, RUNALARM=alarm, ...)`.
   - Creates the solar orbit (MAPL2's `MAPL_Get(MAPL, orbit=orbit, ...)`
     has no MAPL3 replacement - `MAPL_GetOrbit` remains unimplemented,
     see `MAPL_Generic.F90`/`API.F90`). `MAPL_SunOrbitCreateFromConfig`
     (`base/SunOrbit.F90`) still exists and works but takes a legacy
     `type(ESMF_Config)`, unobtainable in MAPL3 - so `Initialize` instead
     reads each orbital parameter directly via `MAPL_GridCompGetResource`
     (labels/defaults copied from `MAPL_SunOrbitCreateFromConfig`'s own
     body, since its `DEFAULT_ORBIT_*`/`DEFAULT_ORB2B_*` constants are
     module-private, not exported) and calls `MAPL_SunOrbitCreate(...)`
     directly with `FIX_SUN=.false.` (no existing resource/precedent for
     this flag anywhere in the repo). The resulting `MAPL_SunOrbit` is
     stored in a new `orbit` field added to `ty_RRTMGP_state` (the
     existing named-private-state type, reused rather than adding a new
     private state) and fetched back in `Run` via
     `_GET_NAMED_PRIVATE_STATE(gc, ty_RRTMGP_state, PRIVATE_STATE, rrtmgp_state)`.
   - `IM`/`JM`/`LM`/`LONS`/`LATS` (previously from `MAPL_Get`): now
     `MAPL_GridCompGet(gc, num_levels=LM, _RC)` +
     `MAPL_GridGet(esmfgrid, IM=IM, JM=JM, _RC)` +
     `MAPL_GridGetCoordinates(esmfgrid, longitudes=LONS, latitudes=LATS, _RC)`,
     matching IRRAD exactly (note: `LONS`/`LATS` changed from `pointer`
     to `allocatable` to match `MAPL_GridGetCoordinates`'s intent(out)
     arrays).
   - `INTERNAL_ESMF_STATE=internal` -> `MAPL_GridCompGetInternalState(gc, internal, _RC)`.
   - All ~34 live `MAPL_TimerOn(MAPL, "NAME", ...)`/`MAPL_TimerOff(MAPL, "NAME", ...)`
     call sites in `Run` -> `MAPL_GridCompTimerStart(gc, "NAME", ...)`/
     `MAPL_GridCompTimerStop(gc, "NAME", ...)` (mechanical, via a
     comment-aware `sed` so the ~14 already-dead/commented-out
     `! call MAPL_TimerOn(MAPL,...)` lines inside the RRTMGP helpers
     were left untouched).
   - The module-scope RRTMGP helpers hoisted in step 2
     (`compute_gas_optics`, `compute_cloud_optics_mcica`,
     `compute_sprlyr_diags_predelta`, `compute_delta_scale`,
     `compute_sprlyr_diags_postdelta`, `compute_rte_sw`,
     `PROCESS_RRTMGP_BLOCK`) each carried a vestigial
     `type(MAPL_MetaComp), intent(inout) :: MAPL` dummy arg used *only*
     by those already-dead/commented-out timer calls (disabled inside
     OMP parallel regions, not thread-safe) - removed the parameter
     entirely from each signature/declaration and from every call site,
     matching IRRAD's equivalent hoisted helpers (which carry no
     `gc`/`MAPL` parameter at all). `shrtwave` was the one exception:
     its `MAPL_TimerOn`/`Off` calls are live (not disabled), so its
     `MAPL` dummy argument was renamed to `gc` (`type(ESMF_GridComp),
     intent(inout) :: gc`) instead of removed, and its sole call site
     (inside `SORADCORE`, host-associated with `gc`) updated to pass
     `gc`.
   - `type(MAPL_MetaComp), pointer :: MAPL` in `Run` and the
     `MAPL_GetObjectFromGC(gc, MAPL, _RC)` call were removed/commented
     out, matching IRRAD's `! type(MAPL_MetaComp), pointer :: MAPL`.
   - **Deferred to step 12, NOT fixed here**: `ImportSpec`/`ExportSpec`/
     `InternalSpec` (`type(MAPL_VarSpec), pointer :: ...(:)`) and the
     `MAPL_VarSpecGet` calls inside the load-balancing block (all
     `MAPL_VarSpec` usage in `Run` is confined to between the
     `"-BALANCE"` timer start/stop) are **confirmed to have zero MAPL3
     equivalent** (`MAPL_VarSpec`/`VarSpecGet` - zero matches anywhere
     in `src/Shared/@MAPL`). Left untouched (still referencing the
     nonexistent type) with a `TODO(step 12)` comment at the
     declarations - this whole block needs a redesign as part of step
     12's load-balancing API investigation, and was already
     non-compilable before this step's changes.

11. **DONE - Edge (`VLOC=E`) 0-based bounds remap.**
   - **Prerequisite fix found and done first**: Solar's `Run` still
     called the old MAPL2 `MAPL_GetPointer(state, ptr, 'NAME', _RC)`
     API at ~140 call sites - confirmed **completely absent from
     MAPL3** (`grep -rl "MAPL_GetPointer\b"` across `src/Shared/@MAPL`
     only turns up the legacy MAPL2-era `apps/mapl_acg.pl` Perl script,
     not any real Fortran symbol or macro). Considered fully switching
     these to the ACG-generated `#include "Solar_DeclarePointer___.h"`/
     `"Solar_GetPointer___.h"` (like IRRAD does), but rejected that:
     unlike IRRAD, most of Solar's manual fetches use a **different**
     local variable name than the field's own short_name (e.g.
     `call MAPL_GetPointer(import, PLL, 'PLE', _RC)`, `RRI` for `'RI'`,
     `ALBIMP` reused across 4 different fields in sequential blocks) -
     switching to ACG's generated declarations (which use the
     short_name as the variable name) would mean renaming every
     downstream numerics usage of `PLL`/`RRI`/etc. throughout ~3000
     lines of `SORADCORE`/`UPDATE_EXPORT`, unverifiable without a build
     (Solar still isn't wired into CMake). Instead did a mechanical
     `sed` rename `MAPL_GetPointer(` -> `MAPL_StateGetPointer(` (real
     MAPL3 function, `superstructure/state/StateGetPointer.F90` -
     confirmed identical positional signature `(state, farrayPtr,
     itemName, unusable, isPresent, rc)`, and confirmed no Solar call
     site uses the old `ALLOC=` keyword, which the new function
     doesn't support) across all ~140 call sites - zero behavior
     change, just the correct MAPL3 name for the same operation. This
     matches how IRRAD itself calls `MAPL_StateGetPointer` directly for
     its own dynamic/per-band fields not covered by its ACG include.
   - **CORRECTION to this step's original text**: the claim that "the
     ACG generator's `emit_declare_pointer()` already emits `contiguous`
     unconditionally" is **wrong** - checked the live
     `apps/MAPL_GridCompSpecs_ACG.py` `emit_declare_pointer()` (and a
     real generated `Irrad_DeclarePointer___.h` in the `ifx/Debug`
     build dir) - ACG emits plain `real(kind=...), pointer ::
     NAME(:,:,:)`, no `contiguous` anywhere. IRRAD's actual remap
     pattern instead declares its own **local** `real, pointer,
     contiguous, dimension(:,:,:) :: p3d` scratch variable in `Run`
     (not ACG-generated) and remaps through it (`p3d => X; X(bounds) =>
     p3d`) - this is what was replicated for Solar.
   - Identified every genuine (non-load-balancing) Edge fetch site by
     checking actual 0-based-index usage, not just the `.rc`'s `VLOC=E`
     tag alone - **`PREF`** (`z`/`E` import, 1D reference pressure) is
     tagged Edge in `Solar_StateSpecs.rc` but is used exclusively with
     plain 1-based indices (`PREF(1)`, `PREF(LM)`, `PREF(K)` for `K` in
     `1..LM`) everywhere in the file, so it was correctly **left
     un-remapped** - remapping it to `0:LM` would have silently shifted
     every access by one level and been a real bug, not a fix. Likewise
     the `case ('PLE')`/`'FSWN'`/etc. references inside `SORADCORE`'s
     load-balancing pack/unpack logic (already `TODO(step 12)`-flagged,
     confirmed broken/deferred) reassign `PLE`/`FSW`/etc. to a
     completely different repacked 2D "daytime-only column" array, not
     the gridded 3D state pointer - correctly left untouched.
   - Added one local `real, pointer, contiguous, dimension(:,:,:) ::
     p3d` scratch declaration in `Run` and a second, separate one in
     `UPDATE_EXPORT` (each contained procedure needs its own scratch
     variable in scope).
   - Remapped, right after each fetch: `AS_PTR_PLE` (`Run`'s AERO/aerosol-
     optics block, 2 fetch sites) and, in `UPDATE_EXPORT`: `PLL` (aliased
     from `'PLE'`), the 8 INTERNAL fields `FSWN`/`FSCN`/`FSWUN`/`FSCUN`/
     `FSWNAN`/`FSCNAN`/`FSWUNAN`/`FSCUNAN` (mandatory/no-COND, no
     `associated()` guard needed, matching IRRAD's INTERNAL-side
     treatment), and the 12 EXPORT fields `FSW`/`FSC`/`FSWNA`/`FSCNA`/
     `FSWD`/`FSCD`/`FSWDNA`/`FSCDNA`/`FSWU`/`FSCU`/`FSWUNA`/`FSCUNA`
     (each guarded by `if (associated(...))`, matching IRRAD's
     EXPORT-side `Update_Flx` treatment, since exports can be
     unassociated if not requested downstream - confirmed Solar's own
     code already assumes this, e.g. `if (associated(FSW)) FSW(:,:,L) =
     ...`).


12. **IN PROGRESS (investigated, not yet implemented) - The
    load-balancing block (`MAPL_LoadBalance`/`MAPL_BalanceWork`) in
    `Run` is unique to Solar (IRRAD has no analog).**
    - **Confirmed working, no changes needed**: `MAPL_BalanceCreate`/
      `MAPL_BalanceWork`/`MAPL_BalanceDestroy` (`mp_utils/
      MAPL_LoadBalance.F90`) still exist in MAPL3 with signatures
      identical to what Solar already calls (checked every real call
      site in `Run` against the actual subroutine signatures). Same for
      the `MAPL_DimsHorzVert`/`MAPL_DimsHorzOnly`/`MAPL_DimsVertOnly`
      dims-classification constants (`mapl_Constants`, reachable via
      plain `use MAPL`).
    - **Small real gap found**: `MAPL_Distribute`/`MAPL_Retrieve` (the
      `Direction=` values Solar passes to `MAPL_BalanceWork`) are
      listed in `MAPL_Public_API.md` but `mp_utils/API.F90`'s `use
      mapl_LoadBalance_mod, only:` list only re-exports the 4 `Balance*`
      functions, not these two parameters - an apparent oversight in
      MAPL3's own aggregator module. Workaround: add a direct `use
      mapl_LoadBalance_mod, only: MAPL_Distribute, MAPL_Retrieve` in
      Solar (bypassing the incomplete umbrella re-export for just these
      two names).
    - **The real blocker**: the pack/unpack machinery (~1000 lines
      across two `"-BALANCE"`-timed regions) generically iterates every
      import/internal field via `MAPL_VarSpecGet(ImportSpec(K),
      DIMS=..., SHORT_NAME=..., _RC)` - and `MAPL_VarSpec`/
      `MAPL_VarSpecGet` are confirmed **completely absent** from MAPL3
      (see step 10's note - re-confirmed here). Checked for a
      metadata-attribute-based MAPL3 replacement (a `'DIMS'`
      `ESMF_Info` attribute automatically attached to fields by
      `MAPL_GridCompAddSpec`, mirroring how `base/NCIO.F90` reads/
      writes a `'DIMS'` info attribute for its own restart/history I/O
      purposes) - confirmed this does NOT exist; `NCIO.F90`'s
      `ESMF_InfoSet(...,'DIMS',...)` calls are internal to that file's
      own I/O logic, not something `gridcomp_add_spec` attaches
      generically to every field it creates. So there is no drop-in
      runtime-introspection replacement for what `MAPL_VarSpecGet` used
      to provide.
    - **Recommended fix (not yet implemented) - CHOSEN APPROACH**:
      instead of a hardcoded name/DIMS table (an earlier idea, now
      superseded), reconstruct the per-field metadata from the ESMF
      states + fields, since MAPL3 does carry all of it on the fields
      themselves (verified against `src/Shared/@MAPL` source). This
      keeps the loop data-driven (no static table to drift out of sync
      with `Solar_StateSpecs.rc`) while dropping `MAPL_VarSpec`
      entirely:
      - **Data source**: drop the `ImportSpec`/`ExportSpec`/
        `InternalSpec` arrays. Pull field lists from the `import`/
        `internal` `ESMF_State` via
        `ESMF_StateGet(state, itemNameList=names)` then
        `ESMF_StateGet(state, name, field)`. Names come free from the
        item list.
      - **`SHORT_NAME`** -> the `itemNameList` entry (or
        `MAPL_FieldGet(field, short_name=)`).
      - **`DIMS`** -> `ESMF_FieldGet(field, rank=)` for slice counting
        (`rank==2` -> 1 slice; `rank==3` -> 3rd extent, the old
        `size(ptr3,3)`).
      - **Vertical-only detection** (`z`/`PREF`, which shares `rank`
        with `xy`) -> `horizontal_dims_spec == HORIZONTAL_DIMS_NONE`
        (from `MAPL_FieldGet(field, horizontal_dims_spec=)`; equiv.
        ESMF `geomDimCount==0`). Confirmed in `MAPL_Generic.F90`
        `gridcomp_add_spec`: `.rc` `DIMS=z` maps to
        `HORIZONTAL_DIMS_NONE`, `xy`/`xyz` -> `HORIZONTAL_DIMS_GEOM`.
      - **Full per-dim shape in one call** (optional convenience):
        `MAPL_FieldGetLocalElementCount(field, local_count, _RC)`
        (public via `infrastructure/esmf/API.F90`) returns the whole
        shape array; no single `ESMF_FieldGet` arg returns the full
        shape.
      - **`UNGRIDDED_DIMS`** count ->
        `MAPL_FieldGet(field, ungridded_dims=ugd)` then
        `ugd%get_num_ungridded()` (`UngriddedDims`,
        `infrastructure/esmf/UngriddedDims.F90`). Handles the
        `ungrd_num_bands_solar` fields (`FSWBANDN`, `DRBANDN`, ...).
      - **`default=def`** -> rely on MAPL3 applying `fill_value` at
        allocation (`FieldClassAspect` calls
        `FieldSet(payload, fill_value)`); drop the `def` arg in unpack
        and verify against a baseline run. Fallback: read `/_FillValue`
        via `ESMF_InfoGet` (`KEY_FILL_VALUE`,
        `utils/MAPL_ESMF_InfoKeys.F90`). For SOLAR this is effectively
        a no-op: only 4 of 159 internals set `FILL` (`TAULOPAR`,
        `TAUMDPAR`, `TAUHIPAR`, `TAUTTPAR`), all to `MAPL_UNDEF`, which
        is already the buffer initialization value.
        - **What `def` actually does** (`INT_VARS_3` unpack loop at
          ~L3746-3800, `UnPackIt` at ~L6860): the load balancer only
          computes **daytime** columns (`daytime`/`MSK` true). When
          unpacking the repacked buffer back into the full gridded
          internal array, daytime cells receive the computed value;
          **nighttime** cells have no computed value, so for
          Out-only internals (`.not. IntInOut(K)`) `UnPackIt`
          overwrites them with `def` (the MAPL2 spec `DEFAULT=`) so
          the whole array is consistent (e.g. zero flux, `MAPL_UNDEF`
          optical depth). For `InOut` internals `def` is deliberately
          NOT passed (`UnPackIt`'s `default` is `optional`; see the
          comment above the loop) so nighttime cells retain their
          previous "aged" values. `MAPL_VarSpecGet(..., default=def)`
          is only called in the non-InOut branch, so `def` is never
          read stale. Hence the replacement only needs to serve the
          Out-only branch: a `fill_value`-based `def`, or dropping
          the arg and trusting allocation-time pre-fill.
      - **Note on `FILL` vs `DEFAULT`**: the step-3 note above lists
        `DEFAULT` as unsupported/dropped - that's the MAPL2 *spelling*.
        The MAPL3 `FILL`/`FILL_VALUE` column DOES work end-to-end
        (`fill_value` is a real `gridcomp_add_spec` dummy arg in
        `MAPL_Generic.F90`), and is what serves the old `default=` role.

      Mapping summary:

      | MAPL2 `MAPL_VarSpecGet` | MAPL3 replacement |
      | --- | --- |
      | `SHORT_NAME` | `ESMF_StateGet(itemNameList=)` |
      | `DIMS` | `ESMF_FieldGet(rank=)` + `horizontal_dims_spec == HORIZONTAL_DIMS_NONE` |
      | `UNGRIDDED_DIMS` | `ungridded_dims%get_num_ungridded()` |
      | `default` | MAPL3 `fill_value` pre-fill (drop arg); fallback `/_FillValue` |

      **Concrete edit plan (spec arrays -> ESMF_State).** Verified there
      are exactly THREE `MAPL_VarSpecGet` calls (L1807, L2033, L3754;
      L769 is only a comment) and that the spec arrays are used only at:
      declarations L770-772, `size()` at L1773-1774, and those three
      `MAPL_VarSpecGet` calls. `ExportSpec` is DECLARED BUT NEVER READ ->
      just delete it.

      1. Declarations (replace L770-772):
         ```fortran
         character(len=ESMF_MAXSTR), allocatable :: ImportNames(:), InternalNames(:)
         type(ESMF_Field) :: field
         ```
         (`internal` ESMF_State local already exists at L762; `import`
         is the Run arg. No `ExportSpec` replacement needed.)

      2. Populate names + counts (near where the spec arrays used to be
         fetched, before the `NumImp = size(...)` block at L1773):
         ```fortran
         call ESMF_StateGet(import, itemCount=NumImp, _RC)
         call ESMF_StateGet(internal, itemCount=NumInt, _RC)
         allocate(ImportNames(NumImp), InternalNames(NumInt), _STAT)
         call ESMF_StateGet(import, itemNameList=ImportNames, _RC)
         call ESMF_StateGet(internal, itemNameList=InternalNames, _RC)
         ```
         Then `size(ImportSpec)` -> `NumImp`, `size(InternalSpec)` ->
         `NumInt` (L1773-1774 become redundant / drop).

      3. Site 1 - input loop (replace the L1807 call; needs DIMS +
         SHORT_NAME):
         ```fortran
         NamesInp(K) = ImportNames(K)
         call ESMF_StateGet(import, ImportNames(K), field, _RC)
         call SolarFieldGetDims(field, DIMS, ugdims, SlicesInp(K), _RC)
         ```
         (`SlicesInp(K)` from the helper replaces the later
         `ESMFL_StateGetPointerToData` + `size(ptr3,3)` slice count for
         the non-aerosol branch; keep the AERO special-case as-is.)

      4. Site 2 - internal loop (replace the L2033 call; needs
         SHORT_NAME + DIMS + UNGRIDDED_DIMS):
         ```fortran
         NamesInt(K) = InternalNames(K)
         call ESMF_StateGet(internal, InternalNames(K), field, _RC)
         call SolarFieldGetDims(field, DIMS, ugDim(K), SlicesInt(K), _RC)
         ```

      5. Site 3 - unpack loop (replace the L3754 `default=def` call):
         ```fortran
         call ESMF_StateGet(internal, InternalNames(K), field, _RC)
         call SolarFieldGetFill(field, def, _RC)
         ```
         (Or drop `def` entirely and rely on MAPL3 `fill_value`
         pre-fill; for SOLAR it's a no-op as noted above.)

      Full replacement mapping:

      | MAPL2 | MAPL3 |
      | --- | --- |
      | `type(MAPL_VarSpec), pointer :: ImportSpec(:)` | `character(len=ESMF_MAXSTR), allocatable :: ImportNames(:)` + `import` state |
      | `type(MAPL_VarSpec), pointer :: InternalSpec(:)` | `character(len=ESMF_MAXSTR), allocatable :: InternalNames(:)` + `internal` state |
      | `type(MAPL_VarSpec), pointer :: ExportSpec(:)` | **delete** (declared, never read) |
      | `size(ImportSpec)` / `size(InternalSpec)` | `ESMF_StateGet(state, itemCount=)` |
      | `MAPL_VarSpecGet(spec, SHORT_NAME=)` | the `itemNameList` entry |
      | `MAPL_VarSpecGet(spec, DIMS=, UNGRIDDED_DIMS=)` | `SolarFieldGetDims(field, ...)` |
      | `MAPL_VarSpecGet(spec, default=)` | `SolarFieldGetFill(field, ...)` or drop |

      **Local helper routine (recommended structure)**: wrap the per-field
      metadata reads in one module-scope private subroutine so the two
      load-balancing loops stay close to their current `select case (DIMS)`
      shape. Place it near the hoisted RRTMGP helpers (after
      `end subroutine Run`); it takes `field` explicitly, no host
      association needed.
      ```fortran
      subroutine SolarFieldGetDims(field, dims, num_ungridded, num_slices, rc)
         type(ESMF_Field), intent(inout) :: field
         integer, intent(out) :: dims          ! MAPL_Dims{HorzVert,HorzOnly,VertOnly}
         integer, intent(out) :: num_ungridded ! old ugDim(K)
         integer, intent(out) :: num_slices    ! 2D slices
         integer, optional, intent(out) :: rc
         integer :: status, rank
         type(HorizontalDimsSpec) :: hspec
         type(MAPL_UngriddedDims) :: ugd
         integer, allocatable :: local_count(:)
         call ESMF_FieldGet(field, rank=rank, _RC)
         call MAPL_FieldGet(field, horizontal_dims_spec=hspec, ungridded_dims=ugd, _RC)
         num_ungridded = ugd%get_num_ungridded()
         if (hspec == HORIZONTAL_DIMS_NONE) then
            dims = MAPL_DimsVertOnly; num_slices = 0 ! z (PREF)
         else if (rank >= 3) then
            dims = MAPL_DimsHorzVert ! xyz
            call MAPL_FieldGetLocalElementCount(field, local_count, _RC)
            num_slices = local_count(3) ! == old size(ptr3,3)
         else
            dims = MAPL_DimsHorzOnly; num_slices = 1 ! xy
         end if
         _RETURN(_SUCCESS)
      end subroutine SolarFieldGetDims
      ```
      Optional companion for the fill-value fallback (only the internal
      unpack loop needs it):
      ```fortran
      subroutine SolarFieldGetFill(field, def, rc)
         type(ESMF_Field), intent(in) :: field
         real, intent(out) :: def
         integer, optional, intent(out) :: rc
         integer :: status
         type(ESMF_Info) :: info
         def = MAPL_UNDEF
         call ESMF_InfoGetFromHost(field, info, _RC)
         if (ESMF_InfoIsPresent(info, KEY_FILL_VALUE, _RC)) &
              call ESMF_InfoGet(info, KEY_FILL_VALUE, def, _RC)
         _RETURN(_SUCCESS)
      end subroutine SolarFieldGetFill
      ```
      Loop usage collapses to:
      ```fortran
      call ESMF_StateGet(import, names(K), field, _RC)
      NamesInp(K) = names(K)
      call SolarFieldGetDims(field, DIMS, ugdims, SlicesInp(K), _RC)
      if (DIMS == MAPL_DimsVertOnly) cycle ! skip PREF
      ```

      **`use`-reachability: all symbols come through the existing
      `use MAPL`** (verified against `src/Shared/@MAPL` API aggregators) -
      no extra `use mapl_*_mod` lines needed:
      - `MAPL_FieldGet` <- `mapl_field_api` (infrastructure/field/API.F90)
      - `MAPL_FieldGetLocalElementCount`, `HorizontalDimsSpec`,
        `HORIZONTAL_DIMS_NONE`, `operator(==)` <- `mapl_esmf_api`
        (infrastructure/esmf/API.F90; the last three via a bare
        `use mapl_HorizontalDimsSpec_mod` there)
      - `UngriddedDims` <- `mapl_esmf_api`, **aliased as
        `MAPL_UngriddedDims`** (use that spelling for the type)
      - `KEY_FILL_VALUE` <- `mapl_utils_api` (utils/API.F90 bare
        `use mapl_esmf_info_keys_mod`)

      One `.rc` cross-check when wiring in: ungridded-only fields
      (`FSWBANDN` = `xy` + `ungrd_num_bands_solar`) come back as
      `MAPL_DimsHorzOnly` with `num_slices=1`; their extra axis is still
      driven by the loop's existing `num_ungridded`/`ugDim` handling,
      matching the current structure.

      Deferred for now at the user's request - pick this up as the next
      concrete task when resuming step 12.

13. **DONE - `Irrad_SetServices`-style external wrapper.** Added a
    standalone `Solar_SetServices(gc, rc)` subroutine after `end
    module`, delegating to the module's `SetServices`, matching
    `Irrad_SetServices` exactly (`use ESMF` + `use
    GEOS_SolarGridCompMod, only: mySetServices => SetServices`, then
    `call mySetServices(gc, rc=rc)`).

14. **DONE - Wire into the parent container** (`GEOS_RadiationGridComp.F90`):
    - Uncommented `use GEOS_SolarGridCompMod, only: solarSetServices =>
      SetServices` and added `call MAPL_GridCompAddChild(gc, "SOLAR",
      solarSetServices, "solar.yaml", _RC)`, mirroring IRRAD's own
      `"irrad.yaml"` literal-filename pattern exactly (no need for the
      inline-`ESMF_HConfigCreate` alternative - IRRAD's precedent shows
      the plain filename overload works fine).
    - Added `MAPL_GridCompAddConnection(gc, src_comp="SOLAR",
      src_names="FSW, FSC, FSWNA, FSCNA", dst_comp="<self>", _RC)` -
      these 4 are the only SOLAR exports the SW/combined formulas below
      actually need (verified against the pre-port `git show 9c9a00f`
      formulas - the last commit before later GEOS-MLT/ML-radiation
      additions that are unrelated new features, not part of this
      port).
    - Added matching `FSW`/`FSC`/`FSWNA`/`FSCNA` IMPORT rows (all
      `VLOC=E`) to `../Radiation_StateSpecs.rc`, mirroring the existing
      IRRAD `FLX`/`FLC`/`FLXA`/`FLA` rows exactly (same `SKIP` restart
      mode, same "connected internally, not for an external coupler"
      comment style). Verified via the real ACG script that the new
      rows parse correctly.
    - Filled in `Run`: declared `FSW`/`FSWCLR`/`FSWNA`/`FSCNA` (aliased
      from `'FSW'`/`'FSC'`/`'FSWNA'`/`'FSCNA'`) alongside the existing
      `FLW`/`FLWCLR`/`FLWNA`/`FLA`, fetched them via
      `MAPL_StateGetPointer`, and remapped them through the same `p3d`
      Edge-remap scratch pointer already used for the LW fields (all 4
      are mandatory/no-COND, forced-allocated by the connection, so no
      `associated()` guard needed - same reasoning as the LW fields).
      Extended the existing LW-only `DMI`-based `if` block to also
      compute `RADSW`/`RADSWC`/`RADSWNA`/`RADSWCNA` (same `DMI`, just
      swapping in the `FSW*` pointers), and added the standalone
      `RADSRF = FSW(:,:,LM) + FLW(:,:,LM)` and
      `DTDT = ((FLW(:,:,0:LM-1)-FLW(:,:,1:LM)) + (FSW(:,:,0:LM-1)-FSW(:,:,1:LM)))
      * (MAPL_GRAV/MAPL_CP)` lines - all formulas taken verbatim from
      the pre-port `git show 9c9a00f:GEOS_RadiationGridComp.F90` (the
      last pre-GEOS-MLT commit, to avoid pulling in the unrelated later
      ML-radiation-blending feature).
    - Re-exported the old `CHILD_ID=SOL` promoted exports (`DRPAR`,
      `DFPAR`, `DRNIR`, `DFNIR`, `DRUVR`, `DFUVR`, `DRPARN`, `DFPARN`,
      `DRNIRN`, `DFNIRN`, `DRUVRN`, `DFUVRN`, `FCLD`, `TAUCLI`,
      `TAUCLW`, `CLDTT`, `ALBEDO`, `FSWBAND`, `FSWBANDNA`) via
      `MAPL_GridCompReexport(gc, src_comp="SOLAR", src_name="...",
      _RC)` (confirmed via `MAPL_Generic.F90`'s `gridcomp_reexport`/
      `ComponentSpec%reexport` and the real `GEOS_SuperdynGridComp.F90`
      precedent that this call **self-registers** the export spec - no
      corresponding `.rc` row is needed, unlike the connected-import
      fields above). `DROBIO`/`DFOBIO` are `COND=SOLAR_TO_OBIO` on
      Solar's side, so guarded the same two reexports behind a fresh
      `MAPL_GridCompGetResource(gc, "USE_OCEANOBIOGEOCHEM", DO_OBIO,
      default=0, _RC)` check in the parent, mirroring Solar's own
      `SetServices` gating logic - reexporting an export that doesn't
      exist on the child side would fail.

15. **DONE - Uncommented `GEOSsolar_GridComp`** in the parent
    `CMakeLists.txt`'s `alldirs` list. `SUBCOMPONENTS`/`DEPENDENCIES`
    already just pass through `alldirs`/`MAPL GEOS_Shared ESMF::ESMF`
    with no explicit per-child entries (matching IRRAD's own precedent
    - it isn't separately listed in `DEPENDENCIES` either), so no other
    change was needed there.

16. **Remaining style cleanup** (lower priority, do after functional
    correctness - step 1 already covers the bulk of macro/indentation
    style): `ESMF_Attribute*`->`ESMF_Info*` (Solar doesn't appear to use
    Attribute get/set based on what's visible, but verify).
