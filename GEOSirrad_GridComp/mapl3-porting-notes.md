# MAPL3 (GEOSgcm/mapl3) conventions learned

- MAPL source IS in this workspace: src/Shared/@MAPL (macros in
  src/Shared/@MAPL/include/MAPL_ErrLog.h, MAPL_private_state.h).
- `__RC__`, `VERIFY_`, `RETURN_`, `ASSERT_` are OLD-style macros but still
  valid/functional in MAPL3 (map to `_RC`/`_VERIFY`/`_RETURN`/`_ASSERT`
  when ANSI_CPP is not defined, which is the case for this build).
  So they are NOT strictly "MAPL2-only" - but GWD/new components use the
  `_RC`/`_STAT`/`_RETURN`/`_ASSERT` style exclusively for consistency.
- `_VERIFY`/`_RETURN`/`_ASSERT` do NOT require an `Iam`/`COMP_NAME`
  traceback variable (that's only needed by the ANSI_CPP macro path).
- `_SET_NAMED_PRIVATE_STATE`/`_GET_NAMED_PRIVATE_STATE` macros
  (MAPL_private_state.h) use `ESMF_InternalStateAdd/Get` - NOT compatible
  with manual `ESMF_UserCompSetInternalState/GetInternalState` pattern.
  Don't mix the two for the same internal-state object.
- ACG (mapl_acg) StateSpecs.rc schema: columns are read by header name,
  not fixed position. Confirmed valid columns beyond the basics:
  `COND` (works in IMPORT/EXPORT/INTERNAL - conditional add, references
  local vars or self% private-state members in scope at #include point),
  `ADD2EXPORT` (in INTERNAL - auto-mirrors an internal field as an export
  of the same name), `UNGRIDDED` (can reference self%-based runtime dims,
  e.g. `self%nbins`). Different components use different column subsets
  (GWD: NAME|ALIAS|UNITS|DIMS|VLOC|RESTART|LONG NAME import; GOCART/
  DynCore drop ALIAS if unused). Drop unused columns rather than leaving
  them blank, matching GOCART/DynCore style.
- mapl_acg() CMake call: `mapl_acg(${this} <Comp>_StateSpecs.rc IMPORT_SPECS EXPORT_SPECS INTERNAL_SPECS [GET_POINTERS DECLARE_POINTERS])`
  generates `<Comp>_Import___.h`, `<Comp>_Export___.h`, `<Comp>_Internal___.h`,
  and optionally `<Comp>_GetPointer___.h`/`<Comp>_DeclarePointer___.h` for Run.
- Genuinely dynamic/runtime-named fields (e.g. IRRAD's RATS_DIAGNOSTICS
  per-tracer exports) can't be expressed in YAML - keep as manual
  MAPL_Add*Spec calls in SetServices.
- fprettify (pip) is a Fortran auto-formatter; installed via a venv here
  since system pip is externally-managed (brew). `-i 3` sets 3-space
  indent, but continuation lines are aligned to the opening paren column,
  NOT a fixed offset - it can't produce a literal "5-space fixed
  continuation" style. For that literal style, must hand-format.
- IRRAD component (GEOS_IrradGridComp.F90) is being ported to MAPL3;
  file is large (~4500 lines) and gets reset to the original MAPL2
  version between sessions/turns (user working in a separate checkout?).
  Re-verify current file state before continuing modernization work.
- IMPORTANT: user asked to save memory updates continuously during this
  porting process (not just at end). After each meaningful edit/decision,
  append a note here (or to session memory for in-progress detail) rather
  than waiting until the end of the conversation.

## ESMF_HConfig (replaces ESMF_Config) - confirmed API (2026-07-17)
Real source-level evidence found in DynCore_GridCompMod.F90 (not just
build artifacts) plus other .F90 sources (GOCART2G, FV_StateMod,
BasicVerticalGrid, generate_scrip_cube):
- Get the hconfig object from a gridcomp: `call MAPL_GridCompGet(gc,
  hconfig=hconfig, ..., _RC)` (declare `type(ESMF_HConfig) :: hconfig`).
  This replaces `ESMF_GridCompGet(gc, config=cf, _RC)` +
  `type(ESMF_Config) :: cf`.
- Scalar existence check: `has_x = ESMF_HConfigIsDefined(hconfig,
  keyString='KEY', _RC)` - must be used as an ASSIGNMENT statement (not
  embedded in an `if (...)` condition), because `_RC`/`__RC__` macros
  expand to two semicolon-separated statements
  (`rc=status);_VERIFY(status)`) which is invalid inside an expression/
  if-condition. Pattern: assign to a logical first, branch on it after.
- Scalar value reads (no `default=` keyword supported by these
  functions - use IsDefined + if/else instead):
  `ESMF_HConfigAsI4(hconfig, keyString='KEY', rc=status)` (integer)
  `ESMF_HConfigAsR4(hconfig, keyString='KEY', rc=status)` (real*4)
  `ESMF_HConfigAsString(hconfig, keyString='KEY', rc=status)` (string)
  `ESMF_HConfigAsLogical(hconfig, keyString='KEY', rc=status)` (logical)
- List/array of strings (replaces the old ESMF_ConfigFindLabel +
  ESMF_ConfigGetLen + loop-of-ESMF_ConfigGetAttribute pattern entirely):
  `xlist = ESMF_HConfigAsStringSeq(hconfig, keyString='KEY',
  stringLen=<len>, _RC)` where `xlist` is declared
  `character(len=<len>), allocatable :: xlist(:)` - assigns the whole
  array in one call (uses Fortran allocatable auto-(re)allocation on
  assignment). Get count via `if (allocated(xlist)) n = size(xlist)`.
  This one function call is a MASSIVE simplification vs the old 4-call
  ESMF_Config loop pattern - always prefer it over manual loops.
- Sub-section access: `sub_hconfig = ESMF_HConfigCreateAt(hconfig,
  keyString='SECTION', rc=status)` then query sub_hconfig the same way.
- Applied in GEOS_IrradGridComp.F90: converted `type(ESMF_Config) :: cf`
  (SetServices) and `type(ESMF_Config) :: CF` (Run) both to
  `type(ESMF_HConfig)`, using MAPL_GridCompGet(...,hconfig=...) instead
  of ESMF_GridCompGet(...,config=...). CO2/CO2_PROVIDER scalar reads use
  IsDefined+AsR4/AsString; RATS_DIAGNOSTICS list reads (both in
  SetServices for export-spec generation AND in Run's LW_Driver
  first-time block) now use a single ESMF_HConfigAsStringSeq call
  instead of the old FindLabel/GetLen/loop pattern.
- The `IAm`/`COMP_NAME` traceback-name pattern (`Iam = "Run"`,
  `Iam = trim(COMP_NAME)//Iam`, etc.) is dead code under the new
  `_RC`/`_VERIFY`/`_RETURN` (and old `VERIFY_`/`RETURN_`) macros since
  ANSI_CPP is not defined - only the ANSI_CPP path uses Iam for
  traceback messages. Removed entirely from IRRAD's `Run`, `LW_Driver`,
  `Update_Flx` (declarations + assignments); kept the underlying
  `ESMF_GridCompGet(gc, grid=..., rc=status)` call (still needed for
  ESMFGRID) but dropped its `name=COMP_NAME` argument since COMP_NAME had
  no other use. 2026-07-17.
- MAPL_GridCompGetResource: confirmed (via GWD .fii build artifact) this
  is exported `use MAPL, only: MAPL_GridCompGet, MAPL_GridCompGetResource`
  from the top-level MAPL module (couldn't find its .F90 definition via
  simple grep - likely template-generated - but usage signature is
  confirmed consistent: `call MAPL_GridCompGetResource(gc, "LABEL", var,
  default=..., rc=status)`). As of 2026-07-17, GEOS_IrradGridComp.F90 has
  ZERO remaining old-style `MAPL_GetResource(...)` calls - all resource
  reads already use `MAPL_GridCompGetResource` (confirmed via grep).
- PREFER `MAPL_GridCompGetResource` over manual
  `ESMF_HConfigIsDefined`+`ESMF_HConfigAsXxx`+if/else when the value is a
  simple scalar with a sensible default - `MAPL_GridCompGetResource(gc,
  "LABEL", var, default=..., _RC)` already handles the "not defined ->
  use default" case internally, so it's a direct one-line replacement for
  that whole IsDefined/if-else block. Applies to CO2 and CO2_PROVIDER
  scalar/string reads in IRRAD's SetServices (simplified 2026-07-17).
  Reserve manual ESMF_HConfigIsDefined+AsXxx only for cases needing
  custom branching logic beyond a simple default value.
- NEW MAPL3 dynamic spec-add API: `MAPL_GridCompAddSpec` (from
  MAPL_Generic, procedure `gridcomp_add_spec`) REPLACES the old
  `MAPL_AddImportSpec`/`MAPL_AddExportSpec`/`MAPL_AddInternalSpec` calls
  for cases needing manual/dynamic spec registration (e.g. runtime-named
  fields that can't go in the static YAML). Confirmed full keyword list
  by reading `gridcomp_add_spec` in
  src/Shared/@MAPL/superstructure/generic/MAPL_Generic.F90:
  `gridcomp` (positional/first arg), `state_intent` (required,
  `ESMF_StateIntent_Flag`: `ESMF_STATEINTENT_IMPORT`/`_EXPORT`/
  `_INTERNAL`), `short_name` (required), then optional: `standard_name`
  (this is what the old `LONG_NAME`/YAML `LONG_NAME` column maps to -
  there is ALSO a separate `long_name` keyword but ACG only emits
  `standard_name`), `units`, `dims` (STRING now, not MAPL_Dims* constant:
  `'xy'` = horz-only, `'xyz'` = horz+vert), `vertical_stagger`
  (`type(VerticalStaggerLoc)`: `MAPL_VERTICAL_STAGGER_NONE`/`_CENTER`/
  `_EDGE` - replaces `VLOCATION=MAPL_VLocationNone/Center/Edge`),
  `ungridded_dim_array` (array of `type(UngriddedDim)`, replaces
  `UNGRIDDED_DIMS=(/n/)` - build each element via
  `MAPL_UngriddedDim(n, name='...', units='...')` then pass
  `ungridded_dim_array=[that_var]`), `add_to_export` (logical, replaces
  old `add2export`/`ADD2EXPORT`), `export_name`, `restart_mode`,
  `fill_value`, `itemtype`, `typekind`, `geom`, `vertical_grid`, etc.
  Real source examples (not just ACG output): GWD's ACG-generated
  `GWD_Internal___.h`/`GWD_Export___.h`, and hand-written calls in
  `src/Shared/@MAPL/gridcomps/orbit/MAPL_OrbGridCompMod.F90` (dynamic
  per-instrument export loop - closest analog to IRRAD's per-tracer RATS
  loop) and `src/Shared/@MAPL/gridcomps/statistics/*.F90`.
  Applied to IRRAD's dynamic RATS-diagnostics export/internal-spec loop
  in SetServices (2026-07-17): converted all `MAPL_AddExportSpec`/
  `MAPL_AddInternalSpec` calls (per-tracer loop + CO2_FIXED/DELT +
  FLXU_RAT/FLXD_RAT/FLX_RAT/DFDTS_RAT/SFCEM_RAT) to
  `MAPL_GridCompAddSpec(gc, state_intent=..., short_name=..., ...)`.
- String-valued `MAPL_GridCompGetResource` outputs should be declared
  `character(len=:), allocatable`, NOT `character(len=ESMF_MAXSTR)` -
  confirmed convention in GWD (`character(len=:), allocatable ::
  gridname` + `call MAPL_GridCompGetResource(gc, 'AGCM.GRIDNAME',
  gridname, _RC)`) and GOCART2G/DU2G (`character(len=:), allocatable ::
  point_emissions_srcfilen` + `default='/dev/null'`) and
  FV_StateMod/CapGridComp. A string literal `default='...'` works fine
  with a deferred-length allocatable actual argument - the subroutine
  auto-allocates to the right length. Applied to IRRAD's `gen_str` in
  SetServices (was `character(len=ESMF_MAXSTR)`, now `character(len=:),
  allocatable`) - 2026-07-17. Note: this only applies to the
  MAPL_GridCompGetResource-fed variable; other `gen_str`/string vars in
  Run/LW_Driver used purely for local string-building (not fed by
  MAPL_GridCompGetResource) were left as fixed-length declarations.

## MAPL2 -> MAPL3 GridComp porting workflow (steps, in order)
This is the checklist being followed to port a legacy GEOS gridded
component (reference: GEOS_GwdGridComp.F90 / DynCore_GridCompMod.F90) to
MAPL3. Apply in this order for each component:
1. Reformat header comments (`!BOP`/`!MODULE:`/`!DESCRIPTION:` style,
   indentation) to match GWD/DynCore conventions.
2. **Create the `<Comp>_StateSpecs.rc` YAML file** (ACG input) by
   extracting the component's `MAPL_AddImportSpec`/`AddExportSpec`/
   `AddInternalSpec` calls out of SetServices into three YAML categories
   (`component:`, `schema_version: 2.0.0`, then `IMPORT`/`EXPORT`/
   `INTERNAL` lists). Match column conventions from an existing
   `*_StateSpecs.rc` (GWD/GOCART2G/DynCore) - drop unused columns, use
   `COND` for conditional fields, `ADD2EXPORT` (INTERNAL only) to mirror
   an internal field as an export, `UNGRIDDED` for runtime self%-based
   dims. Genuinely dynamic/runtime-named fields (names/count only known
   at runtime, e.g. IRRAD's per-tracer RATS_DIAGNOSTICS exports) CANNOT
   go in the YAML - leave those as manual `MAPL_Add*Spec` calls in
   SetServices, with a comment explaining why.
3. Add/verify the `mapl_acg(${this} <Comp>_StateSpecs.rc IMPORT_SPECS
   EXPORT_SPECS INTERNAL_SPECS)` call in `CMakeLists.txt`.
4. In SetServices, replace the manual `MAPL_Add*Spec` calls (that were
   moved to YAML) with `#include "<Comp>_Import___.h"` /
   `_Export___.h` / `_Internal___.h`.
5. Convert private-state management from manual
   `ESMF_UserCompSetInternalState`/`GetInternalState` (+ a wrapper type)
   to `_SET_NAMED_PRIVATE_STATE`/`_GET_NAMED_PRIVATE_STATE` macros -
   convert the setter AND all getters together (not compatible to mix).
6. Convert `MAPL_GetResource(MAPL_obj, ...)` calls to
   `MAPL_GridCompGetResource(gc, "LABEL_NO_COLON", var, default=...,
   _RC)` (no `MAPL_MetaComp` object needed, label has no trailing colon).
7. Convert `ESMF_Config`/`ESMF_ConfigGetAttribute`/`ESMF_ConfigFindLabel`/
   `ESMF_ConfigGetLen` usage to `ESMF_HConfig` (see the ESMF_HConfig
   section above for the concrete API).
8. (Optional/risky, scope with user first) 3-space indent / fixed
   continuation-offset reformatting of the whole file - fprettify can't
   produce a fixed continuation offset, so large legacy Run subroutines
   with no test coverage are high-risk to hand-reformat; confirm scope
   before doing this for anything beyond SetServices.

## Reindenting a huge legacy Run subroutine with fprettify (2026-07-17)
Applied to IRRAD's `Run` (~3300 lines incl. contained LW_Driver/Update_Flx/
compute_lw_*/PROCESS_RRTMGP_LW_BLOCK). Workflow that worked, scoped to just
one subroutine (not the whole file, to avoid touching already-correct
SetServices continuation-line style):
1. `pip install fprettify` into a throwaway venv (`python3 -m venv
   /tmp/fprettify-venv && source .../bin/activate && pip install fprettify`).
2. `sed -n 'START,ENDp' file.F90 > snippet.txt` to extract exactly the
   subroutine text (find START/END lines via `grep -n "^   subroutine
   Run(...)"` / `"^end subroutine RUN"`).
3. Wrap the snippet in a throwaway `module scratch_mod ... contains
   <snippet> end module` with the same `use` statements as the real
   module (fprettify needs syntactically-parseable Fortran, not full
   semantic resolution - dummy/incomplete `use` targets are fine).
4. Run `fprettify -i 3 --disable-whitespace -l 500 wrapped.F90` (in-place).
   `--disable-whitespace` is essential - keeps operator/comma spacing
   untouched and ONLY fixes indentation + blank-line collapsing. `-l 500`
   avoids unwanted line-rewrapping at the default 132-col limit.
5. Diff wrapped-before/after to confirm changes are pure whitespace (no
   token/line content changes beyond blank-line collapsing) before trusting
   it - `diff` showing only reindented lines, no added/removed statements.
6. Extract the reformatted subroutine back out (line range shifts slightly
   because fprettify collapses some blank lines) and splice it back into
   the real file between the untouched prefix/suffix, using sed, not
   replace_string_in_file (files this large are impractical for the
   string-replace tool).
7. Keep a `cp file file.bak` before overwriting, in /tmp - this became the
   restore point for a later "undo" request.
Lesson: doing this via fprettify on the WHOLE file also reformats
continuation-line alignment in already-correct code (e.g. SetServices'
fixed 5-space continuation offset gets replaced with paren-aligned
continuations) - always scope fprettify runs to just the subroutine being
touched, never the whole file, unless the user explicitly wants the whole
file's continuation style changed too.

## Undoing a partial in-progress refactor (2026-07-17)
When asked to "undo the work on X" after a prior (possibly only
partially-visible/truncated) turn had refactored X, the file-on-disk is
the source of truth, not the visible conversation transcript - inspect
the current file structure first (grep for
subroutine/end subroutine/contains) rather than assuming the transcript's
last visible state matches disk. In this session, `LW_Driver` had been
converted from a simple contained subroutine
`LW_Driver(IM,JM,LM,LATS,LONS,RC)` (using host association for
MAPL/gc/import/export/clock/INTERNAL/etc.) to a "standalone-style"
signature with ~40 explicit arguments, but was in fact still nested
inside Run's `contains` (comment in the code claimed "standalone module
procedure" but the structure wasn't actually moved to module scope - an
inconsistent half-done state). To undo: located the pre-refactor
`LW_Driver` body + call site in the `/tmp/*.bak` backup taken during the
earlier reindent task, spliced that back over the refactored version
(call site + subroutine body + closing `end subroutine`), then re-ran the
fprettify reindent workflow (above) on the merged file so the restored
block's indentation matched its now-reindented siblings. Always keep
`/tmp` backups during risky multi-step edits on files with no test
coverage - they're the only reliable rollback point when a request to
"undo" arrives later without git history to fall back on (this file is
not actually tracked by git in this checkout).

## LW_Driver standalone-ification attempted again, skipped (2026-07-21)
Asked again to make `LW_Driver` a standalone module-level routine (see
"Undoing a partial in-progress refactor" above for the first failed
attempt). Before touching anything this time, mapped the real scope:
- `LW_Driver` (552-2311) calls `PROCESS_RRTMGP_LW_BLOCK` (2795-3037, RRTMGP
  branch only, ~line 1676), which calls `compute_lw_aer_optics`,
  `compute_lw_cloud_optics_mcica`, `compute_lw_gas_optics`, `compute_lw_rte`
  (2317-2789). All five are currently siblings inside `Run`'s `contains` -
  moving `LW_Driver` alone breaks its ability to call
  `PROCESS_RRTMGP_LW_BLOCK` (an internal procedure can't call a sibling
  host's other internal procedures) - this is almost certainly the root
  cause of the 2026-07-17 half-done state. `Update_Flx` (3040-3504) is
  separate - `LW_Driver` never calls it, doesn't need to move.
- The bigger blocker: `LW_Driver`'s body (plus the nested RRTMGP cluster)
  operates on **~70 pointer variables** that `Run` fetches once via
  `#include "IRRAD_GetPointer___.h"` (every IRRAD import/internal/export -
  `PLE`, `T`, `Q`, `FLX_INT`, `TS`, `CLDLOLW`, etc, read directly from the
  generated file at
  `build/gfortran/Debug/.../IRRAD_DeclarePointer___.h`), plus
  `USE_RRTMGP`/`USE_RRTMG`/`USE_CHOU`, `gc`, `hconfig`, `NB_IRRAD`,
  `band_output`/`any_band_output`, and the RATS-toggle state (`nRATS`,
  `nameRATS`, `first`, `RATNAMES`). This matches (and explains) the prior
  attempt's "~40 explicit arguments" - the true count is even larger.
- Proposed two options to the user: (a) thread all ~70 pointers +
  scalars through `LW_Driver`'s explicit argument list (faithful 1:1, but
  this is the approach that produced the half-done state last time), or
  (b) pass the state objects (`import`, `export`, `internal`, `gc`) plus
  the small scalar set, and have `LW_Driver` do its own
  `MAPL_StateGetPointer`/`#include "IRRAD_GetPointer___.h"`-style fetches
  internally instead of receiving pre-fetched pointers - small fixed
  signature, no 70-argument surgery.
- **User chose to skip this refactor for now** - no code changes made.
  `LW_Driver` remains nested inside `Run`'s `contains`, exactly as before.
  If revisited: option (b) above is the recommended approach; don't
  attempt the giant-explicit-argument-list version again without a very
  systematic identifier-by-identifier cross-check against `Run`'s full
  scope (host-association bugs here are silent at the call site and only
  surface as compile errors deep in the moved body, which is presumably
  how the 2026-07-17 attempt drifted into an inconsistent state).

## ACG DECLARE_POINTERS / GET_POINTERS feature (2026-07-17)
Confirmed by reading `src/Shared/@MAPL/cmake/mapl_acg.cmake` and
`src/Shared/@MAPL/apps/MAPL_GridCompSpecs_ACG.py` (the `-g`/`-d` flags and
`emit_get_pointers`/`emit_declare_pointers` functions) plus the working
example at `src/Shared/@MAPL/apps/tests/acg3/` (ACG3.F90 +
CMakeLists.txt):
- `mapl_acg(target specs.rc IMPORT_SPECS ... GET_POINTERS
  DECLARE_POINTERS)` (bare keywords, no filename, use ACG default names
  `<Component>_GetPointer___.h` / `<Component>_DeclarePointer___.h`)
  generates two extra include files alongside the existing Import/Export/
  Internal spec ones.
- DECLARE_POINTERS emits `real(kind=ESMF_KIND_R4), pointer ::
  NAME(:,:[,:])` for every FIELD-itemtype spec in the given state
  category/categories; GET_POINTERS emits `call
  MAPL_StateGetPointer(<state>, NAME, 'SHORT_NAME', _RC)` (relies on an
  ambient `status` var for the `_RC` macro, and on a variable literally
  named `import`/`export`/`internal` in scope - Fortran is
  case-insensitive so this matches dummy args `import`/`export` and a
  local `type(ESMF_State) :: INTERNAL` fetched via
  `MAPL_Get(...,INTERNAL_ESMF_STATE=INTERNAL)`).
- The generated variable name is the `ALIAS` column value if present,
  else the `SHORT_NAME` - use an `ALIAS` column in the .rc file whenever
  the desired local Fortran variable name must differ from the
  short_name (e.g. to avoid collision between an IMPORT field and an
  INTERNAL field sharing the same short_name, like IRRAD's 'TS').
- Fields with a `COND` column are AUTOMATICALLY wrapped in
  `if (COND) then <get pointer> else nullify(NAME) end if` by the
  generator - so COND-gated fields (conditionally added to the state)
  are already safe to include in a blanket GET_POINTERS call; no need to
  special-case them.
- CAUTION: any spec whose real runtime object is NOT a plain data pointer
  (e.g. IRRAD's `AERO` import, which is fetched via
  `ESMF_StateGet(IMPORT,'AERO',AERO,...)` into a `type(ESMF_State)`, not
  `MAPL_GetPointer`) must NOT be swept into a blanket
  DECLARE_POINTERS/GET_POINTERS for that category - it has no ITEMTYPE
  column set (defaults to FIELD) so ACG would try to declare/fetch it as
  a real pointer, colliding with the manual `type(ESMF_State) :: AERO`
  declaration (compile error). There's no clean per-field opt-out in the
  current schema (ITEMTYPE only maps 'F'/'V', no 'state' option) - so for
  a category containing such a field, it's safer to leave that whole
  category's pointer-fetching manual and only ACG-ify a clean category.
  Applied 2026-07-17: enabled GET_POINTERS/DECLARE_POINTERS for IRRAD (all
  categories technically generate fine since AERO's COND-less blanket
  inclusion just adds one extra unused/harmless pointer decl+fetch in
  practice - but where an IMPORT-category pointer fetch happens in a
  different subroutine scope than the declare, only replace the
  particular manual block being targeted, don't assume you must convert
  every category at once).
- Test the generator directly before trusting a CMake rebuild: `python3
  src/Shared/@MAPL/apps/MAPL_GridCompSpecs_ACG.py Component_StateSpecs.rc
  -i /tmp/i.h -x /tmp/x.h -p /tmp/p.h -g /tmp/g.h -d /tmp/d.h` runs in
  under a second and lets you inspect the exact generated declare/get
  code before wiring it into the CMakeLists.txt `mapl_acg()` call.
- Applied to IRRAD's INTERNAL category (2026-07-17): added `ALIAS` column
  to `IRRAD_StateSpecs.rc`'s INTERNAL rows (FLX->FLX_INT, SFCEM->SFCEM_INT,
  TS->TS_INT, FLXA/FLXAD/FLXAU->*_INT, etc.; DFDTS* and OLRB*/DOLRB* bands
  left unaliased since their local var names already match short_name),
  added `GET_POINTERS DECLARE_POINTERS` to the `mapl_acg()` call in
  CMakeLists.txt, and replaced the manual `real, pointer, dimension(...)
  :: X_INT` declaration block + the manual `call MAPL_GetPointer(INTERNAL,
  ...)` block in `Run` with `#include "IRRAD_DeclarePointer___.h"` /
  `#include "IRRAD_GetPointer___.h"` respectively.

## Session 2026-07-21: SetServices signature fix, whitespace cleanup, Attribute->Info
- **`SetServices` entry-point signature bug**: had
  `type(ESMF_GridComp), intent(inout) :: gc` and
  `integer, optional, intent(out) :: rc`. Fixed to
  `type(ESMF_GridComp) :: gc` (no intent) and `integer, intent(out) :: rc`
  (non-optional), matching the real `I_SetServices` interface
  (`src/Shared/@MAPL/infrastructure/esmf/ESMF_Interfaces.F90`) and the
  working `GEOS_SuperdynGridComp.F90` reference. This wasn't just a style
  nit: `GEOS_RadiationGridComp.F90` (the sibling composite/container
  component, also ported this session) calls `MAPL_GridCompAddChild(gc,
  "IRRAD", irradSetServices, hconfig, _RC)`, and that call path
  (`mapl_UserSetServices_mod::new_ProcSetServices`) strictly checks the
  passed procedure's signature attribute-for-attribute against
  `I_SetServices` - the old `optional`/`intent(inout)` signature would not
  compile there. `Run`'s signature was already correct (fixed in an
  earlier session/turn).
- **Comment indentation**: fixed 85 comment-only lines that sat at column
  0 instead of matching the indentation of the surrounding code (detected
  by comparing each stray comment's indent against the nearest
  preceding/following code line; a handful of ambiguous cases - e.g.
  commented-out body code right before a `case`/`endif` dedent - resolved
  by reading context rather than blanket rule). Verified whitespace-only
  via `diff` on both files with all leading whitespace stripped.
- **Removed pure dash-separator comment lines** (`!----------` style,
  matching `^[[:space:]]*! ?-+[[:space:]]*$`) - 37 lines removed
  throughout the file. These were purely decorative; removing them left
  the surrounding prose/doc comments intact and readable.
- **Continuation-line indentation**: this file's established convention
  (confirmed in `SetServices`) is a **fixed 5-space offset** from the
  starting statement's own indent for every continuation line (`&`
  continuation) of that statement - NOT alignment to the opening
  parenthesis (which is what an automatic formatter like fprettify would
  do instead). Found 66 of 81 continuation groups violating this
  throughout the file and fixed via a small Python script (adjusts only
  leading whitespace of continuation lines; verified whitespace-only via
  the same strip-and-diff technique). **Gotcha**: a naive "line ends with
  `&`" check breaks when a continuation line has an inline comment after
  the `&` on the same physical line (e.g. `.false. , &!  01` in the
  `band_output_supported` array constructor) - the line does NOT
  end with `&` once you look at the raw trailing character, but it DOES
  still continue the statement. Must strip a trailing `!...` comment
  first (tracking `'`/`"` quote state so a `!` inside a string literal
  isn't mistaken for a comment start) before checking
  `endswith('&')`, or whole continuation groups silently get
  under-detected and only partially reindented.
- **`ESMF_AttributeGet` -> `ESMF_InfoGet`**: converted the 6
  `ESMF_AttributeGet` calls on the `AERO` aerosol-provider state (RH/PLE/
  EXT/SSA/ASY field-name lookups plus the `implements_aerosol_optics_method`
  logical) to the modern Info API, matching the established pattern found
  in `GEOS_OgcmGridComp.F90` and `MBundle_IncrementMod.F90`: fetch a
  handle once via `call ESMF_InfoGetFromHost(AERO, aero_info, _RC)`
  (declared `type(ESMF_Info) :: aero_info`, lowercase per user preference
  even though sibling `AERO`/`AS_*` vars are uppercase), then
  `call ESMF_InfoGet(aero_info, key='...', value=..., _RC)` per lookup
  (keyword is `key=`, not `name=`). Left the one `ESMF_AttributeSet` call
  (`band_for_aerosol_optics`) as-is since only `AttributeGet` was in
  scope - Attribute and Info APIs read/write the same underlying
  attribute storage in ESMF so mixing them on the same host object is
  fine, but flagged to the user as an inconsistency they may want cleaned
  up later (`ESMF_InfoSet(aero_info, key=..., value=..., _RC)`).
- **`LW_Driver`'s RRTMGP module `use` block (554-566) cleaned up**: removed
  the column-alignment padding before every `only:` (e.g.
  `use mo_rte_kind,                only: wp` -> `use mo_rte_kind, only: wp`),
  and collapsed the two multi-line continuations
  (`mo_cloud_sampling`/`mo_optical_props`, each split across 3 lines via
  `&`) onto one line each - safe because this target compiles with
  `-ffree-line-length-none` (confirmed in
  `GEOSirrad_GridComp.dir/flags.make`). Immediately after, asked to
  re-split the two resulting overlong lines (121 and 145 chars) - but via
  **two separate `use mo_X, only: ...` statements for the same module**,
  not `&` continuation. Fortran allows multiple `use ModuleName, only:`
  statements for the same module in one scope; contributions accumulate.
  This is the preferred style here over `&`-continued `only:` lists when a
  line gets long.
- **Lowercased stray uppercase keywords** (`IF`/`THEN`/`ENDIF`/`WHERE`/
  `ALLOCATE`/`DEALLOCATE`) throughout the file - 33 lines, all in the RATS
  first-call-setup block (~924-937), the `REFF`/`*_R` negative-clamping
  `WHERE` statements (~979-983, 1187-1191, 2002-2012), the per-band
  aerosol array `ALLOCATE`/`DEALLOCATE` (~1027-1029, 2297-2299), and two
  `if(...)  THEN` lines in `Update_Flx` (~3209, 3275). Left OpenMP sentinel
  directives (`!$OMP PARALLEL DO` / `!$OMP END PARALLEL DO`, ~1670/1710)
  and the `CALL_LAST` config-label references in header comments alone -
  neither is a real lowercase-able Fortran keyword token (OMP directive
  text lives inside a `!` comment; `CALL_LAST` is a resource-file label
  name, not the `CALL` keyword). Implementation approach worth reusing:
  a comment/string-aware whole-word replacer (tracks `'`/`"` quote state
  and stops at an unquoted `!`) rather than a blind `sed` keyword swap -
  a naive replace would have also matched keyword substrings inside
  identifiers/comments/strings.
- **Renamed `IRRAD_StateSpecs.rc` -> `Irrad_StateSpecs.rc`** (mixed case,
  matching the same rename applied to `RADIATION_StateSpecs.rc` ->
  `Radiation_StateSpecs.rc` in the sibling container). Updated in lockstep:
  `CMakeLists.txt`'s `mapl_acg()` call, and every
  `#include "IRRAD_*___.h"` in `GEOS_IrradGridComp.F90` (`Import`,
  `Export`, `Internal`, `GetPointer`, `DeclarePointer`) -> `Irrad_*___.h`,
  since those generated filenames derive from the specs-file basename, not
  a separate config. Left the `component: IRRAD` field inside the .rc file
  itself unchanged (that's independent metadata, not tied to the filename).
  Older dated entries above this one in this log still say
  `IRRAD_StateSpecs.rc`/`IRRAD_*___.h` - those describe what was true at
  the time they were written and were intentionally left as historical
  record rather than rewritten.
- **`ESMF_AttributeSet` -> `ESMF_InfoSet`**: converted the remaining
  `band_for_aerosol_optics` attribute write (the one left over from the
  earlier `ESMF_AttributeGet`->`ESMF_InfoGet` pass) to
  `call ESMF_InfoSet(aero_info, key='band_for_aerosol_optics',
  value=(OFFSET+band), _RC)`, reusing the same `aero_info` handle. No
  `ESMF_Attribute*` calls remain anywhere in this file.
- **`__STAT__` -> `_STAT`**: replaced all 127 occurrences file-wide
  (`0` remained of the old form afterward; `__RC__` was already `0` before
  this pass, i.e. already fully on `_RC`). Confirmed the two macros are
  functionally identical before doing this (`_STAT` = `_RC_(stat,status)`
  = `stat=status);_VERIFY(status`, vs `__STAT__` = `STAT=STATUS);
  _VERIFY(STATUS` in `MAPL_Exceptions.h` - same expansion modulo
  Fortran's case-insensitivity) - this is a pure style/naming convention
  change, not a behavior change.
- **Added a standalone `Irrad_SetServices(gc, rc)` outside the module**
  (after `end module GEOS_IrradGridCompMod`), matching a real precedent
  found in this repo (`MOM_GEOS5PlugMod.F90`/`MOM6_GEOSPlug.F90`, both of
  which have an external `subroutine SetServices(gc, rc)` after their
  `end module` that just delegates to the module's own `SetServices`).
  Purpose: a module procedure's linker symbol gets module-name-mangled by
  the compiler, so a plugin/DSO-style loader (`dlsym`-based, or MAPL's
  `DSOSetServices` child-add path) can't find it directly - an external,
  file-scope subroutine has a plain, predictable symbol name instead. User
  asked for the routine to be named `Irrad_SetServices` specifically (not
  the generic `SetServices` used by the MOM precedent), and for the same
  treatment in `GEOS_RadiationGridComp.F90` (added there as
  `Radiation_SetServices`, recorded in that file's own porting-notes.md).
  Body is a one-line delegation:
  `use GEOS_IrradGridCompMod, only: mySetServices => SetServices` then
  `call mySetServices(gc, rc=rc)`.
- **CO2_FIXED/DELT statespecs migration attempted, left incomplete**: was
  investigating whether the two fixed-name `MAPL_GridCompAddSpec` calls
  (`CO2_FIXED`, `DELT`, originally ~lines 300-314) could move into
  `Irrad_StateSpecs.rc` like the earlier CO2_FIXED-adjacent discovery
  work. Found they're NOT simply movable as unconditional static specs:
  both sit inside the same `if (n .ne. 0) then ... end if` block (opens
  ~line 231, closes ~line 467, found by matching indentation of `if (n
  .ne. 0) then` against candidate `end if` lines - a bare-`end if` grep
  without indent-matching is misleading here since there are ~10 unrelated
  `end if`s at other indents in between) as the per-tracer dynamic loop,
  where `n = size(nameRATS)` (from the runtime `RATS_DIAGNOSTICS:` config
  list). So while `CO2_FIXED`/`DELT`'s *names* are fixed (unlike the
  loop's `'dOLR_'//trim(nameRATS(i))`-style names), their *presence* is
  conditional on `n /= 0`, and `n`/`nameRATS` are computed (originally
  ~line 225-228) *after* the point where the ACG-generated
  `#include "Irrad_Export___.h"` runs (~line 216) - so a static YAML
  `COND` column referencing `n` wouldn't have anything to reference yet at
  that point. Doing this properly means moving the `nameRATS`/`n`
  computation earlier (before the ACG includes, right after the
  `CO2_PROVIDER` assert block) and adding a `COND: n .ne. 0` column for
  both fields in the `.rc` file. Started this edit, user rejected it
  mid-way (no reason given other than redirecting to a different task) -
  **no changes were made toward this migration**; `CO2_FIXED`/`DELT`
  remain manual `MAPL_GridCompAddSpec` calls exactly where they were, and
  the `.rc` file's existing comment near `CLDLOLW` (claiming all of
  `CO2_FIXED, DELT` etc. "cannot be statically enumerated") is still
  slightly imprecise for these two specifically (their names aren't
  dynamic, only their presence is) but was left unedited too, matching
  "no changes made". If revisited, the plan above (move
  `nameRATS`/`n` earlier + `COND: n .ne. 0`) is believed correct but
  UNVERIFIED (no compile check was done).
- **`logger%info` instead of raw `write(*,*)`**: converted the
  `NUM_BANDS`/`TOTAL_RAD_BANDS` mismatch warning in `Run` (previously
  guarded by `if (MAPL_am_I_Root()) then ... write(*,*) ... end if`) to
  `logger%info(...)` calls, dropping the explicit root-rank guard - no
  real `logger%info` call site found anywhere in this repo
  (`DynCore_GridCompMod.F90`, `AdvCore_GridCompMod.F90`) wraps it in a
  `MAPL_am_I_Root()` check, implying the logger handles rank filtering
  itself. Wiring: added `use pflogger, only: logger_t => logger` at
  module level, `class(logger_t), pointer :: logger` declared in `Run`
  (NOT `SetServices` - first attempt mistakenly added it to
  `SetServices`'s scope before realizing the `write(*,*)` block being
  converted actually lives in `Run`, a different subroutine with its own
  separate locals; had to revert and redo in the right scope), and fetched
  via the existing `MAPL_GridCompGet(gc, grid=..., hconfig=..., ...)` call
  by adding `logger=logger` to it (no separate call needed). **Format
  string caution**: real usage examples in this repo only show `%f` and
  plain-string (no-placeholder) `logger%info` calls - no confirmed
  logical-value format specifier was found (guessed `%l0` once, then
  deliberately backed out that guess rather than ship an unverified format
  token). Used `merge('T','F', some_logical)` string concatenation instead
  to build the message text manually wherever a boolean needed to appear -
  safe because it only relies on the confirmed plain-string call form.
- **`PROCESS_RRTMGP_LW_BLOCK` cluster hoisted to module level**: unlike
  `LW_Driver` (the ~70-pointer problem, see above - still nested inside
  `Run`), the RRTMGP-block-processing cluster it calls -
  `compute_lw_aer_optics`, `compute_lw_cloud_optics_mcica`,
  `compute_lw_gas_optics`, `compute_lw_rte`, `PROCESS_RRTMGP_LW_BLOCK` -
  turned out to be **fully self-contained**: every identifier each one
  references is either an explicit dummy argument, a local declaration,
  its own local `use` statement, or a call to one of these sibling
  helpers. None of them reach into `Run`'s scope at all (verified by
  reading each body in full, not just skimming). So all five moved to
  module level (siblings of `SetServices`/`Initialize`/`Run`) with **zero
  signature changes** - pure relocation + dedent by 3 spaces, one at a
  time (`compute_lw_aer_optics` first, confirmed clean, then the
  remaining four as one contiguous block since they sit back-to-back in
  the file). `LW_Driver` still calls `PROCESS_RRTMGP_LW_BLOCK` fine
  afterward with no changes on its end - a host's contained procedure
  calling a module-level procedure has always been ordinary Fortran, no
  special wiring needed. `Run`'s own `contains` now holds only
  `LW_Driver` and `Update_Flx`. Verification approach for a reorder this
  large: positional `diff` is useless/misleading here (a moved block of
  ~660 lines makes diff show huge apparent insert/delete chunks even
  though nothing but position changed) - instead compare the
  blank-stripped, leading-whitespace-stripped, **sorted** line multiset
  before/after; identical sorted output confirms no line was lost,
  duplicated, or altered, regardless of where it ended up.
- **`Update_Flx` checked for the same standalone treatment, left as-is**:
  unlike the RRTMGP cluster, `Update_Flx` (still nested in `Run`, dummy
  args just `IM, JM, LM, RC`) directly host-associates from `Run`'s
  scope: the state objects `gc`/`import`/`export` themselves (used in
  `MAPL_StateGetPointer(export/import, ...)` and
  `MAPL_GridCompGetResource(gc, ...)` - none of which are its own dummy
  args), the scheme flags `USE_CHOU`/`USE_RRTMGP`/`USE_RRTMG`, and roughly
  35-40 of `Run`'s ~70 ACG-generated internal-state pointers (`TS_INT`,
  the whole `FLX_INT`/`FLXA_INT`/`FLC_INT`/`FLA_INT` family and their
  up/down variants, `DFDTS`/`DFDTSC`/`DFDTSNA`/`DFDTSCNA`, `SFCEM_INT`).
  Found via cross-referencing every identifier in its body (lines
  2310-2774 at the time) against `Run`'s full declared-symbol set (its
  own locals/args plus every name in the ACG-generated
  `Irrad_DeclarePointer___.h`) - same class of problem as `LW_Driver`
  (skipped earlier for the same reason), just a smaller slice of the same
  pointer set. **User chose to leave it as-is** - no changes made,
  `Update_Flx` remains nested inside `Run`'s `contains` alongside
  `LW_Driver`. If revisited, the same "pass state objects, refetch
  pointers internally" approach recommended for `LW_Driver` applies here
  too.
- **Collapsed column-alignment padding in the 5 module-level helpers'
  `use` statements and declarations** (`compute_lw_aer_optics`,
  `compute_lw_cloud_optics_mcica`, `compute_lw_gas_optics`,
  `compute_lw_rte`, `PROCESS_RRTMGP_LW_BLOCK` - the ones hoisted out of
  `Run` earlier this session): the padded style
  (`real(wp), dimension(:,:,:),   intent(inout) :: urand` with extra
  interior spaces to vertically align `intent(...)`/`only:` across
  neighboring lines) was collapsed to single-space throughout - matches
  the earlier `use`-statement cleanup done for `LW_Driver`'s own RRTMGP
  module-use block. `compute_lw_aer_optics` was found already
  single-spaced when this task started (apparently reformatted
  externally between turns - noticed but not investigated further,
  consistent with this file periodically changing outside this session's
  edits per earlier notes), so only the other four needed the collapse
  (82 lines changed total). Also fixed a stray `integer ::status`
  (missing space before the variable name) noticed in
  `compute_lw_aer_optics` while there. Verified whitespace-only via
  `diff` after normalizing all multi-space runs to one space in both the
  before/after copies - identical output confirms no token was altered,
  only interior spacing.
- **Genuine compile-error fix, not a style issue**: `compute_lw_rte`'s
  `use mo_optical_props, only: ty_optical_props_arry, ty_optical_props_1sc`
  was missing the trailing `l` (`ty_optical_props_1sc` vs the real symbol
  `ty_optical_props_1scl`, which the routine's own `select type (...);
  class is (ty_optical_props_1scl)` block already uses) - a genuine typo,
  not something introduced by any whitespace-only script in this session
  (those never touch individual characters within a token). Fixed by
  restoring the trailing `l`. Root cause/origin unknown - possibly
  pre-existing, possibly from whatever external process periodically
  reformats this file between turns (see earlier notes on that).
- **`Update_Flx` cleanup, in two passes**: (1) collapsed the same kind of
  column-alignment padding as the module-level helpers (`real,
  pointer, dimension(:,:  )   :: TSINST` -> `real, pointer,
  dimension(:,:) :: TSINST` - note the extra spaces *inside* the
  `dimension(...)` parens specifically, used to align the closing `)`
  across neighboring lines, needed a `' +\)' -> ')'` regex pass in
  addition to the usual `' {2,}' -> ' '` collapse, or you're left with a
  single dangling space before the paren instead of none) across its
  whole declaration block (47 lines) - and lowercased `intent(IN
  )`/`intent(OUT)` -> `intent(in)`/`intent(out)` in both `Update_Flx`'s
  own signature and (found while sweeping the file for other uppercase
  `intent(IN`/`intent(OUT` occurrences) `LW_Driver`'s signature too,
  matching the lowercase-keyword convention already applied file-wide
  earlier this session. (2) "Combine like terms": merged declarations
  sharing the exact same type+attributes into one statement - the
  12-variable `real, pointer, dimension(:,:,:) :: FLX, FLXA, FLC, ...`
  export-pointer group and the 21-variable `real, pointer,
  dimension(:,:) :: TSREFF, SFCEM, ...` group (wrapped with `&` since
  it's long) each collapsed from ~12-21 separate one-per-line
  declarations into a single statement; `K, LEV_LOW_MID, LEV_MID_HIGH`
  (plain loop/index integers) merged similarly. Deliberately did NOT
  merge: `N` (carries an inline comment, `!<<>> MSL`, that would be
  ambiguous to attach to a merged line), `PRS_LOW_MID`/`PRS_MID_HIGH`
  (each has a distinct explanatory trailing comment - merging would force
  choosing which comment survives or losing one), `TSINST` and
  `FCLD`/`PREF` (same pointer type/rank as the big export groups, but
  sit under separate `! pointer to import` / implicit-import comment
  headers - merging them into the export groups would blur that
  import-vs-export semantic separation that the blank-line/comment
  structure was already conveying), and `STATUS` (kept alone, matching
  the established convention in every other routine in this file where
  `STATUS`/`RC` get their own declaration line, never bundled with
  ordinary locals). Verified no variable name was lost/duplicated by
  tokenizing the before/after declaration block, stripping type/attribute
  keywords (`real`/`integer`/`pointer`/`dimension`/`allocatable`) and the
  intentionally-recased `IN`/`OUT` tokens, and comparing the remaining
  identifier multisets - exact match.
- **`LW_Driver`'s own ~310-line declaration block cleaned up the same
  way** (lines ~556-866: `use` statements, `intent(in)`/`intent(out)`
  args, and the huge locals section). Two passes: (1) collapsed the
  column-alignment padding throughout (94 lines) and additionally
  lowercased stray ALL-CAPS type/attribute keywords found here
  (`REAL`/`INTEGER`/`ALLOCATABLE`/`DIMENSION`/`TARGET`/`POINTER` used
  in all-caps on a handful of lines, e.g. `REAL, ALLOCATABLE,
  DIMENSION(:,:,:,:), target :: TAUA`) to match the lowercase-keyword
  convention already applied file-wide - this file's declarations hadn't
  needed that treatment before since the smaller routines didn't have
  any all-caps type keywords left. (2) Combined a handful of clearly-safe
  like-type groups: `STATUS, loop_status` (unlike the general "`STATUS`
  stays alone" note above - here they were directly adjacent with no
  distinguishing comment on either, so combining them doesn't contradict
  that convention, it just doesn't apply when there's nothing to lose);
  `TAUA, SSAA, ASYA` (shared a `! 4d: dimensioned (...)` comment above,
  not per-variable); the `TAUA_3d`/`SSAA_3d`/`ASYA_3d`/`CWC_3d`/`REFF_3d`
  3D-pointer group; `AS_PTR_3D, AS_PTR_PLE, AS_PTR_T, AS_PTR_Q`;
  `AS_ARR_RH, AS_ARR_PL`; `AEROSOL_EXT, AEROSOL_SSA, AEROSOL_ASY`; and
  the 10-variable 2D export-pointer group (`CLDPRS, CLDTMP, CLDTTLW,
  CLDHILW, CLDMDLW, CLDLOLW, TSREFF, SFCEM, LWS0, DSFDTS`) - kept
  `TAUIR` out of that last group since it's `dimension(:,:,:)` (3D), not
  `(:,:)` like the rest. Left everything else alone: the many
  single-purpose physics locals each carrying a distinct explanatory
  comment (`TAUCRIT`/`PRS_LOW_MID`/`PRS_MID_HIGH`/`LCLDMH`/`LCLDLM`,
  `iceflglw`/`liqflglw`, the `press_ref_min,ptop` /
  `temp_ref_min,tmin` / `temp_ref_max,tmax` trio which are already
  meaningfully paired two-per-line), and `CO2_3d`/`tmp_3d` (each has its
  own `=> null()` initializer plus an inline `<<>> MSL` tag comment).
  Verified via the same tokenize-and-diff approach as `Update_Flx`
  (stripping type/attribute keywords, comparing identifier multisets
  before/after) - exact match, no variable lost or duplicated.

## Edge (VLOC=E) fields need an explicit 0-based bounds remap (2026-07-22)
**Real, confirmed bug, not a style issue.** Traced why: MAPL3's field
creation (`src/Shared/@MAPL/infrastructure/field/FieldCreate.F90`,
`make_bounds_from_args`) always gives the vertical dimension Fortran
bounds `LU_Bound(1, num_levels)` - i.e. `1:num_levels` - regardless of
whether the field is Center- or Edge-staggered. For an Edge field
(`num_levels = LM+1`), that means `MAPL_StateGetPointer` hands back a
pointer bounded `1:LM+1`, **never** `0:LM`. There is no MAPL3-wide
"Edge means 0-based" convention - confirmed by checking
`FV_StateMod.F90` (FVdycoreCubed, already MAPL3-ported): it explicitly
remaps `ak`/`bk` after fetching them (`ak(0:km) => ak1(1:km+1)`,
comment: "AK and BK are expected to be 0-based"), but its OTHER 3D Edge
field, `PE` ("Pressure at layer edges"), is used with plain 1-based
indexing (`PE(:,:,1)`, `PE(:,:,NPZ+1)`) everywhere, no remap at all.
So it's a per-field, per-component decision, not automatic.

IRRAD's code assumes 0-based throughout without ever remapping (e.g. the
comment `! note: PLE(0:LM) but TLEV(1:LM+1)`, and pervasive `(:,:,0)`/
`(:,:,LM)` indexing in `LW_Driver`/`Update_Flx`). Since this build has
`-fcheck=all` enabled, every one of these accesses would have been a
genuine out-of-bounds runtime error the first time the physics actually
ran - it just hadn't surfaced yet because execution hadn't gotten that
far. Found by cross-referencing `Irrad_StateSpecs.rc`'s `VLOC=E` rows
against actual usage, then confirming the fix is *correct* (not just
silencing a crash) via an independent consistency check: `PREF`'s usage
at ~line 991-1008 (`_ASSERT(... > PREF(1) ...)`, `_ASSERT(... <
PREF(LM) ...)`, a `do while` loop starting at `k=1`) only makes physical
sense if `PREF` is 0-based - deliberately skipping index 0 (trivial
near-zero TOA pressure) and treating `LM` as the last/surface index.
Before the fix, `PREF(LM)` would have been the *second-to-last* element
of a `1:LM+1` array, not the surface - this check independently confirms
the remap is the right fix, not just a bounds-check workaround.

**Fix applied** - `ptr(new_bounds) => ptr` self-remap (Fortran pointer
bounds-remapping with the pointer as its own data-target; legal because
the RHS is fully evaluated - target + current association - before the
LHS descriptor is updated, same as any `a => a`-style pointer
reassignment) immediately after each fetch, at every site an Edge
pointer is fetched:
- **`Run`**, right after `#include "Irrad_GetPointer___.h"` (~line 518):
  all 16 INTERNAL Edge fields (`FLX_INT`, `FLC_INT`, `FLA_INT`,
  `FLXD_INT`, `FLXU_INT`, `FLCD_INT`, `FLCU_INT`, `FLAD_INT`,
  `FLAU_INT`, `FLXA_INT`, `FLXAD_INT`, `FLXAU_INT`, `DFDTS`, `DFDTSC`,
  `DFDTSNA`, `DFDTSCNA`) plus the 2 no-COND IMPORT Edge fields (`PLE`
  3D, `PREF` 1D) - all mandatory/always-associated, so no
  `associated()` guard needed. This covers what `LW_Driver` and
  `Update_Flx` both use via host association from `Run`.
- **`Update_Flx`**, right after its own separate `MAPL_StateGetPointer`
  calls for the 12 EXPORT Edge fields (`FLX`, `FLXA`, `FLXD`, `FLXAD`,
  `FLXU`, `FLXAU`, `FLC`, `FLCD`, `FLCU`, `FLA`, `FLAD`, `FLAU`) -
  **does** need `if (associated(X)) X(...) => X` guards, since export
  pointers can be null if not requested downstream by a coupler/HISTORY.
  (`Update_Flx` fetches these itself rather than using `Run`'s
  ACG-fetched copies of the same short names - the two are different
  variable instances due to Fortran scoping/shadowing; `Run`'s own copies
  of the EXPORT-category pointers are fetched by the blanket
  `GET_POINTERS` call but never actually used in `Run`'s own body, only
  in `Update_Flx`'s locally-shadowed ones.)
- **Both RATS-diagnostics fetch sites** (inside `LW_Driver`'s
  `if (nRATS .gt. 0) then` block, and again inside `Update_Flx`'s RAT
  diagnostic section): `DFDTS_RAT`, `FLX_INT_RAT`, `FLXU_INT_RAT`
  (4D, edge is dim 3 of 4, remapped as
  `X(1:IM,1:JM,0:LM,1:nRATS) => X`), plus `FLXD_INT_RAT` at the
  `LW_Driver` site only (not fetched at all in the `Update_Flx` site).
  `SFCEM_INT_RAT` is 3D with no vertical dimension at all (surface-only
  RAT diagnostic) - does NOT need remapping, confirmed by its usage
  (`SFCEM_INT_RAT(:,:,n)`, no `0`/`LM` edge indexing).
- Verified completeness two ways: cross-checked every `VLOC=E` row in
  `Irrad_StateSpecs.rc` against the remap list, and grepped the whole
  file for `X(:,:,0...)`/`X(:,:,LM...)` index patterns to confirm every
  matched variable name was already on the remap list (none missed).

### Follow-up: rank-remapping needs CONTIGUOUS - fixed at the ACG root, not per-callsite
The self-remap pattern above (`X(new_bounds) => X`) failed to compile:
`Error: Rank remapping target must be rank 1 or simply contiguous`.
Fortran's rank-remapping rule requires the data-target to be rank-1 or
provably "simply contiguous" - a plain `pointer` declaration (no
`CONTIGUOUS` attribute) doesn't qualify even when self-referencing.
First fix attempt was a per-callsite workaround: route every remap
through a small `contiguous`-declared scratch pointer
(`p3d => X; X(bounds) => p3d`), needing separate `p3d`/`p1d`/`p4d`
scratch variables in `Run`, `LW_Driver`, and `Update_Flx`. User asked
"can we update ACG instead?" - much better, since this same
rank-remapping need will recur for any other MAPL3 component wanting to
0-base an Edge field, and fixing it once at the generator is a real
root-cause fix vs. N per-component workarounds.

Fixed in `src/Shared/@MAPL/apps/MAPL_GridCompSpecs_ACG.py`,
`emit_declare_pointer()`: added `, contiguous` to the emitted
`pointer` attribute list unconditionally for every ACG-generated
declaration. Confirmed safe first: `compute_rank()` (same file) only
ever returns rank >= 1 for a spec that reaches this function (its
`base_rank` table is `{'z':1, 'xy':2, 'xyz':3}` - there's no rank-0/
scalar case), and `CONTIGUOUS` is invalid for scalars but always legal
for rank>=1 pointers - so the blanket addition can't break a scalar
declaration. This is also just formalizing an assumption MAPL already
relies on elsewhere (`FV_StateMod.F90`'s `ak(0:km) => ak1(1:km+1)`
remap already assumes the underlying ESMF per-DE field data is
contiguous - it just wasn't previously asserted to the compiler via the
attribute). Verified via direct generator invocation (not a full
rebuild): `python3 MAPL_GridCompSpecs_ACG.py Irrad_StateSpecs.rc -d
/tmp/d.h ...` - confirmed `real(kind=ESMF_KIND_R4), pointer, contiguous
:: FLX_INT(:,:,:)` in the output.

After the ACG fix, first reverted `GEOS_IrradGridComp.F90`'s remaps back
to plain direct self-remaps (`X(bounds) => X`, no scratch pointer),
assuming `CONTIGUOUS` on the ACG-generated declaration would be enough.
**That assumption was wrong** - the exact same
"Rank remapping target must be rank 1 or simply contiguous" error
recurred, even with `X` now declared `pointer, contiguous`. Root cause:
`gfortran-15` rejects a *self*-referencing rank remap specifically -
`CONTIGUOUS` on `X` is necessary but not sufficient when `X` is also its
own data-target in the same statement. Tried the obvious-looking
alternative next - add `TARGET` too - but that's outright illegal
Fortran: `gfortran-15` immediately rejects
`real, pointer, contiguous, target :: X(:,:,:)` with "POINTER attribute
conflicts with TARGET attribute" (confirmed via an isolated one-file
test compile, not just reasoning about it). So `POINTER`+`TARGET`
together is never an option, self-remap or not.

**Final, verified-working fix**: keep the ACG `CONTIGUOUS` addition
(still correct and worth having - it's a real improvement, just not
sufficient alone), and reintroduce the two-variable scratch-pointer
routing in `GEOS_IrradGridComp.F90` (`p3d => X; X(bounds) => p3d`,
`p3d`/`p1d`/`p4d` per rank, matching the `ak`/`ak1` two-variable
pattern in `FV_StateMod.F90` - which was never actually a
self-remap to begin with, just looked similar at a glance). Confirmed
this combination (CONTIGUOUS-declared X, remap routed through a
CONTIGUOUS-declared but *distinct* scratch pointer) compiles AND runs
correctly with an isolated `gfortran-15` test program before reapplying
across `Run`, `LW_Driver`, and `Update_Flx` (all 3 needed their
`p3d`/`p1d`/`p4d` scratch declarations put back). One remaining manual
fix outside ACG's reach: `DFDTS_RAT`/`FLX_INT_RAT`/`FLXU_INT_RAT`/
`FLXD_INT_RAT` are hand-declared in `Run` (dynamic RATS-diagnostics
fields added via manual `MAPL_GridCompAddSpec`, not part of the static
`.rc`/ACG output) - added `contiguous` to that one hand-written
declaration line directly; `Update_Flx` shares these via host
association (no separate local declaration), so this covers both the
`LW_Driver` and `Update_Flx` RATS remap sites.

**Lesson for next time**: don't trust "should work now" reasoning about
pointer-remap contiguity rules without an isolated compile test - the
Fortran standard's rules here are subtle enough (rank-1-source vs.
CONTIGUOUS-source vs. self-vs-distinct-target) that a one-file
`gfortran-15 test.f90` check (seconds) is far cheaper than discovering
a wrong assumption after reverting real code back to a broken state.

The ACG `CONTIGUOUS` addition itself still affects every MAPL3
component using `DECLARE_POINTERS`, not just IRRAD - worth knowing if
something elsewhere ever relied on a generated pointer NOT being
contiguous (unlikely - `CONTIGUOUS` only adds a compiler assertion/
optimization hint, it doesn't restrict what the pointer can be assigned
from, and MAPL_StateGetPointer's underlying ESMF field data is
contiguous per-DE in every realistic case).

## Built as a SHARED library (2026-07-22)
Added `TYPE SHARED` to this component's `esma_add_library()` call in
`CMakeLists.txt`, matching the convention used throughout MAPL3 (e.g.
`src/Shared/@MAPL/mapl/CMakeLists.txt`'s own top-level `MAPL` target,
`superstructure/CMakeLists.txt`, `udunits2f/CMakeLists.txt`). Same
change applied to the sibling `GEOS_RadiationGridComp` container's
`CMakeLists.txt`.

## Committed (2026-07-22)
This session's IRRAD work (the `Irrad_SetServices` standalone routine,
the file/variable renames, the whitespace/keyword-case cleanups, the
`Edge`-field 0-based bounds-remap fix, the `TYPE SHARED` build change,
and this notes file) is committed on `feature/pchakrab/port-to-mapl3`
in this repo (commit `bbc9dc1`, bundled together with the matching
`GEOS_RadiationGridComp.F90` fixes in the same commit since both files
live in this one git repo). Not pushed to `origin` yet.

## LW_Driver refresh alarm restored via a real `Initialize` (2026-07-23)
`Run`'s "refresh" block (the full LW transfer calc, `LW_Driver`) had a
commented-out MAPL2-style gate:
`! call MAPL_Get(MAPL, RunAlarm=alarm, _RC)` plus a commented-out
`ESMF_AlarmIsRinging(alarm, ...)` check, with a note that "no MAPL3
`MAPL_GridCompGet` equivalent exists yet for this" - so `LW_Driver` was
running unconditionally every step instead of on its own `IRRAD_DT`
cadence.

`MAPL_GridCompGet`'s public signature (`superstructure/generic/
MAPL_Generic.F90:310-364`, `gridcomp_get`) confirmed there really is no
`clock`/`RunAlarm` output - that part of the old comment was correct.
But MAPL3 doesn't need one: this component didn't have a custom
`Initialize` at all (`SetServices` only registered `Run`, relying on
the generic default for `Initialize`/`Finalize`). Added one, matching
the already-established idiom used throughout this MAPL3 tree for a
component-owned alarm with no natural persistent-state home (see
`DynCore_GridCompMod.F90`'s `replay_shutoff_alarm`, and the `alarm`
lookups in `GEOS_GigatrajGridComp.F90`, `mom_cap.F90`,
`ice_comp_nuopc.F90`): create the alarm **once**, in `Initialize`,
attached to the real `clock` dummy arg (`ESMF_AlarmCreate(name=
"IRRAD_LW_ALARM", clock=clock, ringInterval=..., sticky=.true.,
_RC)`), then in `Run`, just look it up **by name** off the same clock
each call - `ESMF_ClockGetAlarm(clock, alarmname="IRRAD_LW_ALARM",
alarm=lw_alarm, _RC)`. The clock itself is what persists the alarm
across calls; no module-level/`SAVE` Fortran variable or private-state
slot needed for it (unlike the RATS-toggle `first`/`nameRATS` state a
few lines below, which genuinely has nowhere else to live since this
file still has no other one-time-init hook - the new `Initialize` was
deliberately left to do only the alarm setup, not touch that).

Ring interval comes from the `IRRAD_DT` resource (seconds) via
`MAPL_GridCompGetResource(gc, "IRRAD_DT", irrad_dt, _RC)` - no
`default=` given, so it's mandatory, matching how `NUM_BANDS`/`CO2` are
treated elsewhere in this file. Confirmed `IRRAD_DT` is a real,
already-configured key (not something new invented for this) - it
shows up as `IRRAD_DT: 3600` in the example `AGCM.rc` configs under
`src/Applications/@UMD_Etc/UMD_rc/` and the MIT-plug instructions dir.

`Run`'s refresh block itself is now just the originally-commented-out
code with the comment markers stripped and `alarm` renamed `lw_alarm`
(declared as a plain, non-`SAVE` local - it's re-fetched from the clock
every call, so no `SAVE` is needed here either):
```fortran
if (ESMF_AlarmIsRinging(lw_alarm, _RC)) then
   call ESMF_AlarmRingerOff(lw_alarm, _RC)
   call MAPL_GridCompTimerStart(gc, "LW_DRIVER", _RC)
   call LW_Driver(IM, JM, LM, lats, lons, _RC)
   call MAPL_GridCompTimerStop(gc, "LW_DRIVER", _RC)
endif
```

**Not build-verified** - this session had no way to compile the file in
isolation (the full dependency chain is heavy); build it before
trusting it compiles/links. Also not yet committed (left unstaged
alongside an unrelated, pre-existing `NUM_BANDS` error-message refactor
already sitting in this file's working tree - see git history/status,
not this session's alarm work).
