# GEOS_SolarGridComp MAPL3 porting plan

Status: in progress. `GEOSsolar_GridComp` is still commented out of the
parent container's `alldirs` (see `../CMakeLists.txt`) and
`GEOS_SolarGridComp.F90` is still the unported MAPL2 version. Steps 1
and 2 below are done (see branch `feature/pchakrab/port-solar-to-mapl3`).

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
3. **Extract state specs into `Solar_StateSpecs.rc`.** Pull every
   `MAPL_AddImportSpec`/`AddExportSpec`/`AddInternalSpec` call out of
   `GEOS_SolarGridComp.F90`'s `SetServices` into YAML (`IMPORT`/`EXPORT`/
   `INTERNAL` categories, `ALIAS` column for internal-vs-import name
   collisions). Solar is much larger than Irrad here - it has the huge
   `#ifdef SOLAR_RADVAL` block of internals (TAUxxPAR, COTxxPAR, etc.)
   and the OBIO-conditional exports (`SOLAR_TO_OBIO`), so expect many
   `COND` columns (`USE_RRTMG .or. USE_RRTMGP`, `SOLAR_TO_OBIO`,
   `#ifdef SOLAR_RADVAL` fields stay as manual `MAPL_GridCompAddSpec`
   calls guarded by the same `#ifdef` since ACG can't express
   compile-time conditionals). The `OSRBbbRG`/`ISRBbbRG`/`TBRBbbRG`
   per-band dynamic exports (loop over `ibnd`) are genuinely dynamic -
   leave those as manual `MAPL_GridCompAddSpec` calls, same treatment as
   IRRAD's RATS diagnostics loop.
4. **Add `mapl_acg()` to `CMakeLists.txt`** with `IMPORT_SPECS
   EXPORT_SPECS INTERNAL_SPECS GET_POINTERS DECLARE_POINTERS`, and add
   `TYPE SHARED` to `esma_add_library()`.
5. **Replace static `MAPL_Add*Spec` blocks** in `SetServices` with
   `#include "Solar_Import___.h"` / `_Export___.h` / `_Internal___.h`.
6. **Convert the RRTMGP internal state.** Solar's `ty_RRTMGP_state`/
   `ty_RRTMGP_wrap` uses manual `ESMF_UserCompSetInternalState`/
   `GetInternalState` - convert to `_SET_NAMED_PRIVATE_STATE`/
   `_GET_NAMED_PRIVATE_STATE` macros (don't mix the two patterns).
7. **`MAPL_GetResource(MAPL,...)` -> `MAPL_GridCompGetResource(gc,
   "LABEL", var, default=..., _RC)`** (drop trailing colon on labels)
   throughout `SetServices` and `Run`.
8. **`ESMF_Config` -> `ESMF_HConfig`** wherever Solar reads config
   directly (if any beyond `MAPL_GridCompGetResource`).
9. **Fix entry-point signatures**: `SetServices(gc, rc)`/
   `Run(gc, import, export, clock, rc)` with **no `intent`/`optional`**
   on `gc`, and `rc` non-optional `intent(out)`. This is required
   before `MAPL_GridCompAddChild` can add Solar as a child from
   `GEOS_RadiationGridComp.F90`.
10. **Check `MAPL_MetaComp`/`MAPL_Get`/`MAPL_TimerOn` usage** in `Run` -
   Solar's `Run` has heavy use of `MAPL_Get(MAPL, ...)` for grid/clock
   info and `MAPL_TimerOn(MAPL,"TOTAL")`/`"PRELIMS"` - convert to
   `MAPL_GridCompGet`/`MAPL_GridCompTimerStart`/`Stop`.
11. **Edge (`VLOC=E`) 0-based bounds remap - Solar needs this too,
   likely worse than IRRAD.** Solar's SORAD core assumes `PLE(0:LM)`-
   style indexing pervasively (fluxes at layer interfaces indexed
   top-down from L=0). Every Edge-staggered import/internal/export
   pointer (fluxes `FSW`/`FSC`/etc., `PLE`) needs the same `contiguous`
   scratch-pointer remap pattern IRRAD used (`p3d => X;
   X(1:IM,1:JM,0:LM) => p3d`) - the ACG generator's
   `emit_declare_pointer()` already emits `contiguous` unconditionally
   now (fixed during the IRRAD port), so no generator change needed,
   just apply the remap at each fetch site in `Run`/`SORADCORE`/
   `Update_Flx`.
12. **The load-balancing block (`MAPL_LoadBalance`/`MAPL_BalanceWork`)**
    in `Run` is unique to Solar (IRRAD has no analog) - verify these
    MAPL APIs still exist unchanged in MAPL3's `MAPL_Generic`/
    `MAPL_LoadBalanceMod`; this is new territory not covered by the
    IRRAD port.
13. **`Irrad_SetServices`-style external wrapper.** Add a standalone
    `Solar_SetServices(gc, rc)` subroutine after `end module`,
    delegating to the module's `SetServices`, matching
    `Irrad_SetServices`/`Radiation_SetServices`.
14. **Wire into the parent container** (`GEOS_RadiationGridComp.F90`):
    - Uncomment `use GEOS_SolarGridCompMod, only: solarSetServices =>
      SetServices`.
    - `call MAPL_GridCompAddChild(gc, "SOLAR", solarSetServices,
      "solar.yaml", _RC)` (or inline `ESMF_HConfigCreate(content='{}',
      ...)` if Solar needs no per-child yaml - check whether Solar
      reads anything via its own hconfig vs. resource file).
    - Add `MAPL_GridCompAddConnection` pulling Solar's needed exports
      (`FSW`/`FSC`/etc., whatever `RADSW`/`RADSWC`/`RADSWNA`/
      `RADSWCNA`/`DTDT`/`RADSRF` need) into `<self>`.
    - Add matching rows to `../Radiation_StateSpecs.rc`'s IMPORT
      category for those connected fields (mirroring the existing IRRAD
      `FLX`/`FLC`/... rows).
    - Fill in the SW/combined pointer fetches + `RADSW`/`RADSWC`/
      `RADSWNA`/`RADSWCNA`/`DTDT`/`RADSRF` calculations in `Run`, using
      the pre-MAPL3 formulas (`git log` on the pre-port
      `GEOS_SolarGridComp.F90`/old `GEOS_RadiationGridComp.F90` has
      them - the parent's porting notes explicitly flag this as the
      follow-up work).
    - Re-export the old `CHILD_ID=SOL` promoted exports (`DRPAR`,
      `FCLD`, `ALBEDO`, `TAUCLI`, etc.) via `MAPL_GridCompReexport(gc,
      src_comp="SOLAR", src_name="...", _RC)`.
15. **Uncomment `GEOSsolar_GridComp`** in the parent `CMakeLists.txt`'s
    `alldirs` list, and add `SOLAR` to `SUBCOMPONENTS`/`DEPENDENCIES` if
    it isn't automatically picked up.
16. **Remaining style cleanup** (lower priority, do after functional
    correctness - step 1 already covers the bulk of macro/indentation
    style): `ESMF_Attribute*`->`ESMF_Info*` (Solar doesn't appear to use
    Attribute get/set based on what's visible, but verify).
