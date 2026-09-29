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
  - `DEFAULT` (**CORRECTED**): the MAPL2 column name is unsupported, but
    the functional replacement is the MAPL3 **`FILL`** column ->
    `fill_value=` on `MAPL_GridCompAddSpec`, applied at allocation by
    `FieldClassAspect%allocate`. The original "inert/droppable" reading
    here caused step 3 to drop 9 of Solar's 13 `DEFAULT=MAPL_UNDEF`
    internals - since fixed (see step 12's `default=def` note). Rule:
    every MAPL2 `DEFAULT=x` must become `FILL=x`. Note also that MAPL2
    zero-filled internals with no `DEFAULT=`; MAPL3 leaves them
    uninitialised when `FILL` is empty. To keep bootstrapped starts
    zero-diff with MAPL2 before the first REFRESH, every Solar INTERNAL
    row now has an explicit `FILL` (146 x `0.0`, 13 x `MAPL_UNDEF`);
    see the comment block above the INTERNAL table in
    `Solar_StateSpecs.rc`.
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
     - Re-verified 2026-09-26 against `base/SunOrbit.F90` (develop):
       `MAPL_SunOrbitCreate` takes 19 positional args (`CLOCK,
       ECCENTRICITY, OBLIQUITY, PERIHELION, EQUINOX, EOT, ORBIT_ANAL2B,
       ORB2B_YEARLEN, ORB2B_REF_YYYYMMDD, ORB2B_REF_HHMMSS,
       ORB2B_ECC_REF, ORB2B_ECC_RATE, ORB2B_OBQ_REF, ORB2B_OBQ_RATE,
       ORB2B_LAMBDAP_REF, ORB2B_LAMBDAP_RATE, ORB2B_EQUINOX_YYYYMMDD,
       ORB2B_EQUINOX_HHMMSS`) + optional `FIX_SUN` (defaults `.false.`
       when absent) + `RC`; Solar's call (L697-L705) matches in order
       and type. The `MAPL_GridCompGetResource` defaults (L678-L695)
       equal `DEFAULT_ORBIT_*`/`DEFAULT_ORB2B_*` (0.0167, 23.45, 102.0,
       80, 365.2596, 20000101, 115856, 0.016710, -4.2e-5, 23.44,
       -1.3e-2, 282.947, 1.7195, 20000320, 73500) and `EOT`/
       `ORBIT_ANAL2B` default `.false.` as in
       `MAPL_SunOrbitCreateFromConfig`. Labels are the MAPL2 ones minus
       the trailing `:`.
     - `MAPL_SunOrbit`, `MAPL_SunOrbitCreate`, `MAPL_SunGetInsolation`
       (generic: `SOLAR_1D`/`SOLAR_2D`/`SOLAR_ARR_INT`),
       `MAPL_SunGetSolarConstant`, `MAPL_SunGetLocalSolarHourAngle` are
       all re-exported by `base/API.F90`, so plain `use MAPL` suffices.
     - The two `MAPL_SunGetInsolation` call sites (`SORADCORE` ~L1696
       with `INTV=TINT, currTIME=, TIME=SUNFLAG, DIST=`; `UPDATE_EXPORT`
       ~L4029 with `INTV=DELT, clock=, TIME=SUNFLAG, ZTHN=...`) resolve
       to `SOLAR_2D` and are unchanged from MAPL2 - no port needed.
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


12. **DONE (compiles with ifx) - The
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
      - **`default=def`** -> **CORRECTED 2026-09-25; the two earlier
        ideas below are both WRONG**: (a) "rely on MAPL3 allocation-time
        `fill_value` pre-fill and drop `def`" and (b) "fallback: read
        `/_FillValue` via `ESMF_InfoGet`". Verified against MAPL2
        v2.71.0 (`~/workspace/code/GEOSgcm/main/src/Shared/@MAPL`) and
        MAPL3 source:
        - **MAPL2 zero-filled every internal with no `DEFAULT=`**
          (`generic/VarSpec.F90:256` `usableDEFAULT=0.0`;
          `base/Base/Base_Base_implementation.F90:159` `init_value =
          0.0`). **MAPL3 does NOT**: `FieldClassAspect%allocate`
          (`superstructure/generic/specs/FieldClassAspect.F90:254`)
          only calls `FieldSet(payload, fill_value)` `if
          (allocated(this%fill_value))`; otherwise the payload is left
          uninitialised. So the 155 non-`FILL` internals cannot get
          their `def` from MAPL3 - it must be an explicit `0.0`.
        - **`fill_value` is not runtime-introspectable**: nothing in
          MAPL3 writes `KEY_FILL_VALUE` (`/_FillValue`) onto a field's
          `ESMF_Info` (`FieldInfo.F90` only imports the key; the only
          `_FillValue` writers are NCIO/Mesh file metadata). The
          `ESMF_InfoGet` fallback would always miss. Drop
          `SolarFieldGetFill` from the plan.
        - **Step 3 silently lost 9 defaults**: the MAPL2 source had 13
          `AddInternalSpec`s with `DEFAULT=MAPL_UNDEF` (`COSZSW`,
          `CLDTTSW`, `CLDHISW`, `CLDMDSW`, `CLDLOSW`, `COTLOPAR`,
          `COTMDPAR`, `COTHIPAR`, `COTTTPAR`, `TAULOPAR`, `TAUMDPAR`,
          `TAUHIPAR`, `TAUTTPAR`) but `Solar_StateSpecs.rc` only carried
          `FILL=MAPL_UNDEF` on the 4 `TAU*PAR` rows. **FIXED**: added
          `FILL=MAPL_UNDEF` to the other 9 rows (re-ran the real ACG:
          13 `fill_value=MAPL_UNDEF` emissions, names match).
        - **Chosen replacement for site 3**: a module-level private
          `pure function solar_internal_default(name) result(def)` -
          `select case (trim(name))` returning `MAPL_UNDEF` for the 13
          names above, `0.0` otherwise - called as `def =
          solar_internal_default(NamesInt(K))` in the Out-only branch.
          Keep passing `def` to `UnPackIt` exactly as today. The `.rc`
          `FILL` column stays as the documented source of truth and
          must be kept in sync with the `select case` by hand (add a
          one-line comment at both sites saying so). Rejected
          alternative: `FILL=0.0` on all 155 rows - noisy, and still
          wouldn't make the value queryable.
        - **`MAPL_UNDEF` changed sign**: MAPL2 `MAPL_UNDEFINED_REAL =
          huge(1.)` (`shared/Constants/InternalConstants.F90:20`);
          MAPL3 `= -huge(1.)` (`utils/Constants/InternalConstants.F90:18`).
          Solar has zero hard-coded undef literals (`grep huge|1e15`
          empty) and ACG passes `fill_value=MAPL_UNDEF` textually, so a
          MAPL3 build is self-consistent (`where (RI == MAPL_UNDEF)`,
          `BufOut = MAPL_UNDEF`, `def`, exports all agree). Exposure is
          only data crossing the build boundary: (i) MAPL2-written
          Solar internal restarts / ExtData carrying `+huge` must be
          regenerated from a MAPL3 run (or converted `+huge -> -huge`
          once); (ii) baseline validation against MAPL2 output needs an
          undef-aware comparison (treat `+huge` and `-huge` as equal);
          (iii) never introduce a literal - always the symbolic
          constant.
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
      | `default` | `solar_internal_default(name)` lookup (see CORRECTED note above) |

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
         def = solar_internal_default(NamesInt(K))
         ```
         (No field/state access needed. Do NOT drop `def`: MAPL3 does
         not zero-fill, see the CORRECTED note above.)

      Full replacement mapping:

      | MAPL2 | MAPL3 |
      | --- | --- |
      | `type(MAPL_VarSpec), pointer :: ImportSpec(:)` | `character(len=ESMF_MAXSTR), allocatable :: ImportNames(:)` + `import` state |
      | `type(MAPL_VarSpec), pointer :: InternalSpec(:)` | `character(len=ESMF_MAXSTR), allocatable :: InternalNames(:)` + `internal` state |
      | `type(MAPL_VarSpec), pointer :: ExportSpec(:)` | **delete** (declared, never read) |
      | `size(ImportSpec)` / `size(InternalSpec)` | `ESMF_StateGet(state, itemCount=)` |
      | `MAPL_VarSpecGet(spec, SHORT_NAME=)` | the `itemNameList` entry |
      | `MAPL_VarSpecGet(spec, DIMS=, UNGRIDDED_DIMS=)` | `SolarFieldGetDims(field, ...)` |
      | `MAPL_VarSpecGet(spec, default=)` | `solar_internal_default(NamesInt(K))` |

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
      Companion for the unpack loop's `def` (replaces the withdrawn
      `SolarFieldGetFill`; see CORRECTED note above):
      ```fortran
      ! Keep in sync with the FILL column of Solar_StateSpecs.rc.
      pure function solar_internal_default(name) result(def)
         character(*), intent(in) :: name
         real :: def
         select case (trim(name))
         case ('COSZSW', 'CLDTTSW', 'CLDHISW', 'CLDMDSW', 'CLDLOSW', &
               'COTLOPAR', 'COTMDPAR', 'COTHIPAR', 'COTTTPAR', &
               'TAULOPAR', 'TAUMDPAR', 'TAUHIPAR', 'TAUTTPAR')
            def = MAPL_UNDEF
         case default
            def = 0.0 ! MAPL2's implicit default for internals
         end select
      end function solar_internal_default
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
      - (`KEY_FILL_VALUE` no longer needed - fallback withdrawn.)

      One `.rc` cross-check when wiring in: ungridded-only fields
      (`FSWBANDN` = `xy` + `ungrd_num_bands_solar`) come back as
      `MAPL_DimsHorzOnly` with `num_slices=1`; their extra axis is still
      driven by the loop's existing `num_ungridded`/`ugDim` handling,
      matching the current structure.

      Deferred for now at the user's request - pick this up as the next
      concrete task when resuming step 12.

      **Status: DONE.** `GEOS_SolarGridComp.F90` now builds against
      MAPL3 with ifx (`make GEOSsolar_GridComp` in `build/ifx/Debug`;
      only pre-existing unused-variable remarks). What was changed:
      - `use mapl_LoadBalance_mod, only: MAPL_Distribute, MAPL_Retrieve`
        and `use mapl_HorizontalDimsSpec_mod, only: HorizontalDimsSpec,
        HORIZONTAL_DIMS_NONE, operator(==)` added to the module.
      - `ImportSpec`/`ExportSpec`/`InternalSpec` (`MAPL_VarSpec`)
        replaced by `ImportNames(:)`/`InternalNames(:)` filled via
        `ESMF_StateGet(state, itemCount=)` + `itemNameList=`.
      - New module helper `SolarFieldGetDims(field, dims,
        num_ungridded, ug_extent, rc)` (after `end subroutine Run`)
        rebuilds the VarSpec `DIMS`/`UNGRIDDED_DIMS` view from
        `MAPL_FieldGet(horizontal_dims_spec=, ungridded_dims=,
        num_levels=)`.
      - New module helper `solar_internal_default(name)` (pure) returns
        `MAPL_UNDEF` for the 13 rows with `FILL: MAPL_UNDEF` in
        `Solar_StateSpecs.rc`, `0.0` otherwise; replaces
        `MAPL_VarSpecGet(default=)` in `INT_VARS_3`. Must be kept in
        sync with the `.rc` FILL column.
      - `INPUT_VARS_1` checks `itemType`; `AERO` (nested state) keeps
        `DIMS = MAPL_DimsHorzOnly` and goes down the existing branch.
      - `INT_VARS_1`: `if (associated(ugdims))` -> `if (num_ungridded >
        0)`, `ugDim(K) = ug_extent`. Slice counts still come from
        `size(ptr,3)` on the state pointers (zero-diff).
      - All 13 `ESMFL_StateGetPointerToData(` -> `MAPL_StateGetPointer(`.
      Run-testing in progress via `regression/solar-sa` (see step 18);
      still needs the regression comparison against MAPL2 once the full
      model links.

      Facts found while re-verifying, on top of the plan above:
      - `HorizontalDimsSpec`/`HORIZONTAL_DIMS_NONE`/`operator(==)` are
        NOT reachable via `use MAPL` (esmf/API.F90 `use`s the module
        but never `public ::`s them) - add `use
        mapl_HorizontalDimsSpec_mod, only: HorizontalDimsSpec,
        HORIZONTAL_DIMS_NONE, operator(==)` next to the
        `mapl_LoadBalance_mod` line.
      - `ESMFL_StateGetPointerToData` (13 call sites, L1844-L3790) does
        not exist in MAPL3; drop-in is `MAPL_StateGetPointer(state,
        ptr, name, _RC)` (same arg order).
      - The `import` state holds `AERO` as a nested `ESMF_State`
        (SetServices `itemtype=MAPL_STATEITEM_STATE`), so the
        `itemNameList` walk must check `itemType` and route `AERO`
        to the existing special case instead of the field helper.
      - `MAPL_FieldGet(ungridded_dims=)` EXCLUDES the vertical dim;
        use `num_levels=` (>0 -> `MAPL_DimsHorzVert`) rather than
        `rank` to tell `xyz` from `xy`+ungridded (both rank 3).
        `ugDim(K)` = `ugd%get_ith_dim_spec(1)%get_extent()` when
        `get_num_ungridded()==1`, else 0.

      **Optional follow-on (analysed, deferred): extract the two
      `"-BALANCE"` blocks out of `SORADCORE` into module-level
      routines.** SUPERSEDED by step 17 (which generalises this to the
      whole of `SORADCORE`; the balance extraction is its sub-steps
      17.1-17.3, and step 17 should now run BEFORE this step's
      `MAPL_VarSpec` rewrite). Original findings kept for reference:
      - Extents: distribute block ~L1726-L2598 (~873 lines, from the
        `"-BALANCE"` timer start through the second
        `MAPL_BalanceWork(Direction=MAPL_Distribute)` and `NCOL =
        size(Q,1)`); retrieve block ~L3729-L3811 (~83 lines,
        `MAPL_BalanceWork(Direction=MAPL_Retrieve)` + `INT_VARS_3`
        unpack + `MAPL_BalanceDestroy`).
      - The retrieve block is fully generic and moves cleanly.
      - The distribute block is NOT self-contained: the generic work
        (slice counting, `BufInp`/`BufInOut`/`BufOut` allocation,
        `PackIt`, `MAPL_BalanceWork`) is interleaved with two large
        `select case (NamesInp(K))` / `select case (NamesInt(K))`
        blocks (~L1961-L2016 and ~L2199-L2583) that bind ~200 named
        pointer locals of `SORADCORE` (`PLE`, `T`, `Q`, ..., `FSW`,
        `OSRBRGN(k)%p`, every `SOLAR_RADVAL` diagnostic) to
        rank-remapped slices of the buffers. Those pointers are what
        the Chou/RRTMG/RRTMGP code consumes, so they cannot move to a
        module routine without ~200 pointer dummy args.
      - Host-associated state from `Run` that would have to become
        explicit arguments: `gc`, `import`, `internal`,
        `ImportSpec`/`InternalSpec` (or their MAPL3 replacement),
        `LATS`, `AEROSOL_EXT`/`SSA`/`ASY`, `implements_aerosol_optics`,
        `NUM_BANDS_SOLAR`, `USE_RRTMG`, `USE_RRTMGP`, `band_output`,
        `SolarBalanceHandle`, `DYCORE`, `ibnd`, `bb`.
      - State shared between the two blocks (created in distribute,
        consumed in retrieve): `NumMax`, `Num2do`, `daytime`,
        `HorzDims`, `SlicesInp`/`SlicesInt`, `NamesInp`/`NamesInt`,
        `IntInOut`, `rgDim`, `ugDim`, `BufInp`/`BufInOut`/`BufOut`.
      - Proposed design (zero-diff: identical pack order and buffer
        layout):
        - a module-level `type SolarBalance` holding all the shared
          state above plus per-variable buffer offsets `OffInp(:)`/
          `OffInt(:)`;
        - `solar_balance_distribute(gc, import, internal, <names/
          metadata source>, ZTH, SLR, LATS, Ig, Jg, AEROSOL_*,
          include_aerosols, ..., bal, rc)` - the generic part only;
          records offsets instead of binding pointers;
        - `solar_balance_retrieve(gc, internal, <metadata source>,
          LoadBalance, bal, rc)` - the whole retrieve block;
        - type-bound accessors `bal%inp2d(name)` / `bal%out2d(name)`
          (and 1D/3D variants as needed) returning a rank-remapped
          pointer into the buffer, so the two `select case` blocks in
          `SORADCORE` collapse to one-liners (`PLE => bal%inp2d('PLE')`,
          `FSW => bal%out2d('FSWN')`) and stay in `SORADCORE` where the
          pointers live.
      - Payoff: `SORADCORE` shrinks by ~700 lines, the pointer bindings
        lose all loop/offset bookkeeping, and the MAPL3
        `ESMF_StateGet`-based metadata walk from the plan above only
        has to be written once, inside the two new routines.

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

17. **TODO - Slim `SORADCORE` down to a ~100-line driver by hoisting
    everything else to module scope.** Analysed 2026-09-25, not started.
    `SORADCORE` is ~L1349-L3816 (~2470 lines); ~2100 of those can move.
    This supersedes the "Optional follow-on" sketch under step 12 (the
    load-balancing extraction is sub-steps 17.1-17.3 below) and should
    be done BEFORE the step-12 `MAPL_VarSpec` rewrite, so that rewrite
    lands entirely inside two small routines. Solar is now wired into
    CMake (step 15), so every sub-step can be compile-checked, unlike the
    step-2 hoists. Zero-diff throughout: same pack order, buffer layout,
    and call arguments.

    **Why nothing downstream of the balance block can move today**: all
    three schemes (Chou/RRTMG/RRTMGP) consume ~30 input pointers (`PLE`,
    `T`, `Q`, `QL`, ..., `ALBVR`, `ZT`, `SLR1D`, `Ig1D`, `Jg1D`, `ALAT`)
    and write ~200 output pointers (`FSW`, `FSC`, ..., `OSRBRGN(:)`,
    every `SOLAR_RADVAL` diagnostic), all `SORADCORE` locals bound by
    the two `select case` blocks (~L1961-L2016, ~L2199-L2583). The
    enabler is to move those pointers into a derived type.

    **Supporting types (module scope, private)**:
    - `type SolarColumns` - every input/output column pointer listed
      above, incl. `type(rptr1d_wrap) :: OSRBRGN(nbndsw), ISRBRGN(nbndsw)`
      and the whole `#ifdef SOLAR_RADVAL` pointer block (~100 decl lines
      leave `SORADCORE`). The existing `rrtmg_sw(...)` /
      `PROCESS_RRTMGP_BLOCK(...)` calls keep their signatures; args just
      become `cols%COTDTP, ...`.
    - `type SolarBalance` - load-balancing shared state: `handle`,
      `NumMax`, `Num2do`, `NumLit`, `HorzDims(2)`, `daytime(:,:)`,
      `NamesInp/NamesInt`, `SlicesInp/SlicesInt`, `IntInOut`, `rgDim`,
      `ugDim`, per-variable offsets `OffInp(:)/OffInt(:)`,
      `BufInp/BufInOut/BufOut` (allocatable, target), plus the three
      `BUFIMP_AEROSOL_*` pointers.
    - `type SolarWork` - scheme-independent aux arrays, allocatable
      components: `PL`, `RH`, `PLhPa`, `QQ3`, `RR3`, `O3`, `taua`,
      `ssaa`, `asya`. Auto-deallocated on scope exit, replacing the two
      explicit `deallocate` blocks.
    - `type SolarConfig` - the ~20 `Run`-local scalars both schemes
      read: `LM`, `NCOL`, `CO2`, `DIST`, `DOY`, `MG`, `SB`, `SC`,
      `LCLDLM`, `LCLDMH`, `include_aerosols`,
      `implements_aerosol_optics`, `SOLAR_TO_OBIO`, `USE_RRTMG`,
      `USE_RRTMGP`, `band_output(:)`, `SolCycFileName`, `USE_NRLSSI2`,
      `IM_World`, `num_aero_vars`, `currTIME`, `alarm`. Filled once in
      `SORADCORE`. (Alternative is ~20 explicit args on each scheme
      routine; the type is cleaner. IRRAD's `LW_Driver` didn't need
      this because it has one scheme path.)

    **Sub-steps (do in this order; each is a separate commit)**:

    | # | New module routine | Source lines | Notes |
    | --- | --- | --- | --- |
    | 17.0 | define the four types | - | plus `cols` local in `SORADCORE`; no behaviour change |
    | 17.1 | `solar_bind_columns(bal, cols)` | the two `select case` blocks (~450) | pure `bal` -> `cols` pointer binding; the ONLY sub-step touching pointer targets, do first |
    | 17.2 | `solar_balance_distribute(gc, import, internal, <spec/names>, ZTH, SLR, LATS, Ig, Jg, AEROSOL_*, cfg, bal, rc)` | L1726-L2598 minus 17.1 (~420) | `"-BALANCE"` start through 2nd `MAPL_BalanceWork(Distribute)`; records offsets instead of binding pointers |
    | 17.3 | `solar_balance_retrieve(gc, internal, <spec/names>, LoadBalance, bal, rc)` | L3729-L3811 (~83) | already fully generic |
    | 17.4 | `solar_get_insolation(alarm, orbit, LONS, LATS, currTIME, SUNFLAG, SC, esmfgrid, ZTH, SLR, DIST, Ig, Jg, rc)` | L1690-L1720 (~35) | `MAPL_SunGetInsolation` + `SLR*SC` + `Ig/Jg` fill |
    | 17.5 | `solar_prepare_aux(cols, cfg, bal, work, rc)` | L2602-L2691 (~90) | `RH`/`PL`/`PLhPa`/`QQ3`/`RR3`/`O3`/aerosol copy; reads `bal%BUFIMP_AEROSOL_*` |
    | 17.6 | `run_rrtmgp_sw(gc, cfg, cols, work, rc)` | L2715-L3243 (~530) | incl. `_GET_NAMED_PRIVATE_STATE` (needs only `gc`) and the `TEST_` macro define/undef |
    | 17.7 | `run_rrtmg_sw(gc, cfg, cols, work, rc)` | L3245-L3711 (~467) | |

    Optional later sub-splits once 17.6/17.7 compile: `rrtmgp_setup_inputs`
    (mu0/tsi/albedos/p_lay/t_lay/kluges/dzmid, ~120) and `rrtmgp_post`
    (flux normalisation + output load + super-band sums, ~110);
    `rrtmg_flip_inputs` (the `--RRTMG_FLIP` block, ~150) and
    `rrtmg_unflip_outputs` (unflip + `CLDTS` + `COTTP` `where` blocks, ~80).

    **Target shape of `SORADCORE` after 17.7**:
    ```fortran
    subroutine SORADCORE(IM, JM, LM, include_aerosols, currTIME, MaxPasses, LoadBalance, rc)
       ! ~40 decl lines: ZTH, SLR, Ig, Jg, DIST, bal, cols, work, cfg, NCOL
       call solar_get_insolation(..., ZTH, SLR, DIST, Ig, Jg, _RC)
       call solar_balance_distribute(..., bal, _RC)
       call solar_bind_columns(bal, cols)
       NCOL = size(cols%Q, 1)
       cols%COSZSW = cols%ZT
       if (.not. include_aerosols) then ! alias to the "A" internals
          cols%FSW => cols%FSWA; cols%FSC => cols%FSCA; ...
       end if
       cfg = SolarConfig(...)
       call solar_prepare_aux(cols, cfg, bal, work, _RC)
       if (USE_CHOU) then
          call shrtwave(... cols%... , gc, _RC)
       else if (USE_RRTMGP) then
          call run_rrtmgp_sw(gc, cfg, cols, work, _RC)
       else if (USE_RRTMG) then
          call run_rrtmg_sw(gc, cfg, cols, work, _RC)
       else
          _FAIL('unknown SW radiation scheme!')
       end if
       call solar_balance_retrieve(..., bal, _RC)
       _RETURN(_SUCCESS)
    end subroutine SORADCORE
    ```

    **Caveats**:
    - The `FSW => FSWA` (etc.) aliasing at ~L2610 must stay in
      `SORADCORE` (or at the end of `solar_bind_columns`): it happens
      after binding and before the schemes run.
    - `RADSW_BINARY_CLOUDS` resource read + `where (CL > 0.) CL = 1.`
      (~L2622) belongs in 17.5.
    - Step 12's `MAPL_VarSpec` -> `ESMF_StateGet` rewrite then touches
      only 17.2 and 17.3.

18. **IN PROGRESS - Standalone run-testing via `regression/solar-sa`**
    (`mpirun --n 6 ../install/bin/GEOS.x mapl.yaml` from
    `build/ifx/Debug/solar-sa`; failures show up as a `FAIL at line=`
    traceback in `log.run`, with the Solar line at the top). Fixes so far,
    2026-09-28:
    - `MAPL_GridCompSetEntryPoint` calls in `SetServices` now pass
      `phase_name='initialize'` / `phase_name='run'`, matching the
      convention in `docs/mapl2-to-mapl3-port.md`.
    - **`ESMF_GridCompGet(gc, GRID=esmfgrid)` fails at runtime in MAPL3**
      (first `Run` call, old L857): the user-level `gc` has no ESMF grid
      attached - the geom lives on the outer meta component. Replaced
      with `MAPL_GridCompGet(gc, num_levels=LM, grid=esmfgrid, _RC)`
      (`gridcomp_get` in `MAPL_Generic.F90` accepts optional `geom=` /
      `grid=`), which is what IRRAD's `Run` already does. `ESMF_GridCompGet`
      is still fine for `NAME=` only. Same rule applies to any other
      ported component: never ask ESMF for the grid, ask MAPL.

    Fixes 2026-09-29 (run now completes; results match MAPL2 baseline,
    see "Comparison status" below):
    - **Run alarm never rang, so REFRESH (`SORADCORE`) never executed.**
      `solar_run_alarm` was created with only `ringInterval`, and ESMF
      then first rings at `currTime + interval`. MAPL2's generic main
      alarm (`handle_clock_and_main_alarm`, `MAPL_Generic.F90`
      L1328-1411) did much more, and MAPL3 provides nothing equivalent
      for component alarms, so `Initialize` now reproduces it verbatim:
      `RUN_AT_INTERVAL_START` (default `.false.`), `REFERENCE_DATE`
      (default yyyymmdd of currTime) / `REFERENCE_TIME` (default 0) ->
      `ring_time`; if `ring_time > currTime` step it back by whole
      intervals; unless `run_at_interval_start`, back off one clock dt
      (clock advances AFTER Run); walk forward `while ring_time <
      currTime`; `ESMF_AlarmCreate(..., ringTime=, ringInterval=,
      sticky=.false.)`; `ESMF_AlarmRingerOn` if `ring_time == currTime`.
      `Run` keeps its `ESMF_AlarmRingerOff` as in MAPL2. IRRAD has the
      same latent bug - see the TODO in its porting notes.
    - **`ESMF_Attribute*` on the AERO state fails (status 57).** MAPL3
      providers (and `FakeGOCART`) set state attributes with `ESMF_Info`,
      so `Run`'s AERO block now does `ESMF_InfoGetFromHost(AERO,
      aero_info)` once, then `ESMF_InfoGet(aero_info, key=, value=)` for
      `implements_aerosol_optics_method`, `*_for_aerosol_optics`,
      `extinction_in_air_due_to_ambient_aerosol`, etc., and
      `ESMF_InfoSet(aero_info, key='band_for_aerosol_optics', ...)`.
    - **`MAPL_GridGetInterior` no longer exists** (link error).
      Replaced in `SORADCORE` with `MAPL_GridGet(esmfgrid,
      interior=interior, _RC)` (allocatable `integer :: interior(:)`,
      `iBeg/iEnd/jBeg/jEnd = interior(1:4)`), as IRRAD already does.
    - **Stat 151 (already allocated) at the second `SORADCORE` call.**
      `ImportNames`/`InternalNames` are `Run`-scope allocatables that
      `SORADCORE` allocates each call, and it is called twice per REFRESH
      when `do_no_aero_calc`. Added `deallocate(ImportNames,
      InternalNames, _STAT)` to `SORADCORE`'s cleanup.
    - **Edge field `PREF` is 1-based in MAPL3, 0-based in MAPL2**, and the
      `LCLDMH`/`LCLDLM` level search starts at `K=1`, so both super-layer
      boundaries were off by one (122/142 instead of 121/141) and
      `COT{NUM,DEN}{HI,MD,LO}PAR` / `CLD{HI,MD,LO}SW` were wrong. Same
      remap IRRAD uses: `p1d => PREF; PREF(0:LM) => p1d` right after the
      `MAPL_StateGetPointer`. Check every `VLOC=E` pointer this way.

    Run procedure (from `build/ifx/Debug/solar-sa`):
    - Rebuild/install: `cd build/ifx/Debug && bash -c 'module load
      ifx-stack && make -j8 install'`.
    - Each run advances `cap_restart.yaml` to 22:20; reset first with
      `printf 'currTime: 2000-04-14T22:00:00\nrepeatCount: 0\n' >
      cap_restart.yaml`.
    - `bash -c 'module load ifx-stack; export
      LD_LIBRARY_PATH=$PWD/../install/lib:$LD_LIBRARY_PATH; mpirun --n 6
      ../install/bin/GEOS.x mapl.yaml > run.log 2>&1'`. Without the
      `LD_LIBRARY_PATH` export the shared `libGEOSradiation_GridComp.so`
      is "not found" (FAIL in `UserSetServices.F90`).
    - Checkpoints land in `checkpoints/2000-04-14T22:20:00/` (`last`
      symlink); compare with `cmpchk.py` (needs `module load ifx-stack`
      for netCDF4). The script treats a cell as undef if the baseline is
      `_FillValue`-masked or either side is `|huge|` (MAPL2 `MAPL_UNDEF`
      is `+huge`, MAPL3 is `-huge`) and reports undef-mask mismatches
      separately.

    Comparison status vs `/home/pchakrab/input/solar/C12L181/after/
    new-style/SOLAR_{import,internal,export}_after_runPhase1.nc`:
    **bit-for-bit identical in all three states** (2026-09-29), once the
    MAPL2 baseline was captured with a configuration matching the
    standalone:
    - `AERO_PROVIDER: none` in `AGCM.rc` (default `GOCART2G`; options
      `GOCART2G, MAM, none`). `GEOS_ChemGridComp` then creates an empty
      `AERO` state with `implements_aerosol_optics_method = .false.`,
      which is what `FakeGOCART` does here. With the original GOCART2G
      baseline, `FSWN/FSCN/FSWUN/FSCUN/DR*/DF*/FSWBANDN` differed
      O(1e-5 - 7e-3) on day cells (current `FSWN == FSWNAN`).
    - `HISTORY.rc` requesting `[IO]SRB{08..11}RG` from `SOLAR` (only
      bands 08-11 are in `band_output_supported`). The standalone's
      `activate_all_exports: true` turns `band_output` on for those
      bands, so `[IO]SRBbbRGN` are computed here; in the original
      baseline they were never written and held 0 or load-balance
      buffer `MAPL_UNDEF` garbage on day cells.
    - After re-capturing, **also refresh the internal restart** in
      `checkpoints/2000-04-14T22:00:00/SOLAR_internal.nc` from the new
      `before/new-style/SOLAR_internal_before_runPhase1.nc`. With
      `CALLED_LAST` defaulting to 1, `UPDATE_EXPORT` runs before the
      refresh and `[IO]SRBbbRG = [IO]SRBbbRGN * SLR` comes from the
      restart; a stale restart gave `1e15 * SLR` values and `-huge`
      (all-zero faces) in the `*RG` exports.
    - The remaining `RI`/`RL` "undef mismatch" (baseline `_FillValue`
      masked outside cloud, current literal `1e15`) is a file-encoding
      artifact that `cmpchk.py` now treats as undef on both sides.
