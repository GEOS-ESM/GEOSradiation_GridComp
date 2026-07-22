# GEOS_RadiationGridComp MAPL3 porting notes

Status as of this session: `GEOS_RadiationGridComp.F90` (the composite
container that used to own SOLAR + IRRAD + SATSIM children) has been ported
to MAPL3, wired to its **only currently-ported child, IRRAD**. SOLAR and
SATSIM remain commented out / unbuilt (`CMakeLists.txt`'s `alldirs` already
had them commented as "NOT ported to MAPL3 yet" before this session).

## Files touched
- `Radiation_StateSpecs.rc` (new; originally created as
  `RADIATION_StateSpecs.rc`, renamed for consistency) - ACG input for this container's own
  Import/Export specs.
- `CMakeLists.txt` - added the `mapl_acg()` call for the new specs file.
- `GEOS_RadiationGridComp.F90` - full rewrite (see below).

## Scope decision: LW-only for now
Solar isn't wired, so any export that needs both LW (IRRAD) and SW (Solar)
flux can't be correctly computed yet:
- **Populated in `Run`**: `RADLW`, `RADLWC`, `RADLWNA`, `RADLWCNA`, `ALW`,
  `BLW` - these only need IRRAD's exports.
- **Declared in the export spec but left unpopulated**: `DTDT`, `RADSRF`
  (need FLW **and** FSW), `RADSW`, `RADSWC`, `RADSWNA`, `RADSWCNA` (need FSW
  only, which doesn't exist at all without Solar). When Solar is ported,
  add the `MAPL_GridCompAddChild`/`MAPL_GridCompAddConnection` calls for it
  (mirroring the IRRAD pattern below) and fill in the SW/combined pointer
  fetches + calculations in `Run`, matching the original MAPL2 file's logic
  (`git log` on the pre-port version has the full formulas).
- The Solar/Satsim `CHILD_ID=SOL` promoted exports from the old file
  (`DRPAR`, `FCLD`, `ALBEDO`, `TAUCLI`, etc.) are not in the new state spec
  at all yet - add them back via `MAPL_GridCompReexport(gc, src_comp="SOLAR",
  src_name="...", _RC)` once Solar exists as a child.

## Key MAPL3 patterns/discoveries (verified against MAPL source, not guessed)

- **Entry-point subroutine signatures must match `I_SetServices`/`I_Run`
  exactly**: `type(ESMF_GridComp) :: gc` etc. with **no** `intent()`
  qualifier, and `integer, intent(out) :: rc` - **not** `optional`. Found
  via `src/Shared/@MAPL/infrastructure/esmf/ESMF_Interfaces.F90`. Confirmed
  by comparison with the real, working `GEOS_SuperdynGridComp.F90`. (Note:
  `GEOS_IrradGridComp.F90`'s `Run` had `optional` on `rc` and it silently
  compiled anyway when passed to `MAPL_GridCompSetEntryPoint` - but it
  reliably **fails to compile** the moment you try to pass a component's
  `SetServices` as a child via `MAPL_GridCompAddChild`, because that path
  goes through an explicit-interface wrapper
  (`mapl_UserSetServices_mod::new_ProcSetServices`) that strictly checks
  attribute-for-attribute against `I_SetServices`. Always declare `rc`
  non-optional on `SetServices`/`Initialize`/`Run`.)
- **Adding a child with no child-specific config**: `MAPL_GridCompAddChild`
  always needs an `hconfig` (or `hconfig_file`) argument in its public
  overloads. Every real ported GridComp in this repo supplies a `.yaml`
  file (e.g. SuperdynGridComp's `dyn.yaml`). For a child that needs no
  config of its own, an inline empty hconfig avoids adding a placeholder
  file: `hconfig = ESMF_HConfigCreate(content='{}', _RC)` then
  `call MAPL_GridCompAddChild(gc, "IRRAD", irradSetServices, hconfig, _RC)`
  then `call ESMF_HConfigDestroy(hconfig, _RC)`. (There IS a lower-level
  hconfig-free path via `mapl_ChildSpec_mod`/`mapl_UserSetServices_mod`,
  but those are used-but-not-re-exported in
  `src/Shared/@MAPL/superstructure/generic/API.F90` - i.e. deliberately not
  part of the public `use MAPL` surface. Don't reach around the umbrella
  module for this; use the inline-hconfig approach instead.)
- **Parent Run reading a child's export state** (old MAPL2 pattern:
  `MAPL_Get(MAPL, GEX=GEX, ...)` + `MAPL_GetPointer(GEX(child), ..., 
  alloc=.TRUE.)`) has **no direct MAPL3 replacement** - nothing in this
  repo does it. The correct MAPL3-idiomatic replacement is a **self-target
  connection**: declare the needed field(s) as one of the container's own
  IMPORTs in the StateSpecs.rc, then in `SetServices`:
  ```fortran
  call MAPL_GridCompAddConnection(gc, &
       src_comp="IRRAD", src_names="FLX, FLC, FLXA, FLA, DSFDTS0, SFCEM0, TSREFF", &
       dst_comp="<self>", _RC)
  ```
  `"<self>"` is a real sentinel (`StateRegistry.F90`'s `SELF = "<self>"`,
  also used internally by `MAPL_GridCompReexport`'s `reexport()` in the
  opposite direction: child-export -> **my own export**). It is not
  exported as a named Fortran constant anywhere - use the literal string.
  After the connection, `Run` just does a normal
  `MAPL_StateGetPointer(import, ptr, 'FLX', _RC)` on its own `import`
  state; no `alloc=` keyword exists in MAPL3's `MAPL_StateGetPointer` (a
  field simply comes back unassociated if not COMPLETE - but a declared
  connection forces allocation, so this is safe for connected fields).
- **Running a single child in Run**: `MAPL_GridCompRunChild(gc, "IRRAD",
  _RC)` - `phase_name` is optional here and defaults to `'GENERIC::RUN_USER'`
  (the phase a child gets when its own `SetServices` registers `Run`
  without an explicit `phase_name=`, as IRRAD does). The plural
  `MAPL_GridCompRunChildren(gc, phase_name=..., _RC)` requires an explicit
  `phase_name` in its public wrapper even though the underlying method
  allows it to be optional - prefer the singular call per child unless you
  genuinely want "run every child with this exact phase name".
- **`MAPL_MetaComp`/`MAPL_GetObjectFromGC`/`MAPL_Get(MAPL,...)`/
  `MAPL_TimerOn`/`MAPL_GetResource(MAPL,...)`** are all gone from this file,
  replaced by `MAPL_GridCompGet`/`MAPL_GridCompGetResource`/
  `MAPL_GridCompTimerStart`/`MAPL_GridCompTimerStop` taking `gc` directly.
  `MAPL_GridCompGetResource` labels have **no trailing colon** (old style
  was `"LABEL:"`, new is `"LABEL"`).
- `MAPL_GenericSetServices`/`MAPL_GenericInitialize`/
  `MAPL_GenericRunChildren` are not called at all in MAPL3 - the
  `OuterMetaComponent` wrapper handles generic init/child recursion
  implicitly. A custom `Initialize` entry point is still registered here
  (via `MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_INITIALIZE, Initialize,
  _RC)`) only because this component has real one-time setup work
  (`set_inhomogeneity`, `initialize_cloud_subcol_gen`) beyond what the
  generic layer does automatically.

## Edge (VLOC=E) import pointers needed the same 0-based remap as IRRAD (2026-07-22)
Same bug class as the one found and fixed in `GEOS_IrradGridComp.F90`
(see that file's own `mapl3-porting-notes.md` for the full root-cause
writeup): MAPL always creates Edge-staggered fields with Fortran bounds
`1:LM+1`, never `0:LM`, but `Run`'s arithmetic assumes 0-based throughout
(`PLE(:,:,1:LM)-PLE(:,:,0:LM-1)`, `FLW(:,:,0:LM-1) - FLW(:,:,1:LM)`,
etc). Checked `Radiation_StateSpecs.rc`: only the 5 IMPORT fields are
`VLOC=E` - `PLEINST` (this container's own import) and `FLX`/`FLC`/
`FLXA`/`FLA` (connected in from IRRAD). All 9 EXPORT fields (`DTDT`,
`RADLW`, ..., `ALW`, `BLW`) are `VLOC=C` or `VLOC=N` - no export-side
remap needed here, unlike IRRAD which had 12 Edge exports.

These 5 pointers are fetched via plain, hand-written
`MAPL_StateGetPointer` calls in `Run` (this file doesn't use ACG's
`GET_POINTERS`/`DECLARE_POINTERS` for `Run` at all, only the `_Import___.h`/
`_Export___.h` spec-registration includes in `SetServices`), so the ACG
`CONTIGUOUS` fix doesn't reach them automatically - added
`contiguous` directly to their hand-written declarations (`PLE`, `FLW`,
`FLWCLR`, `FLWNA`, `FLA`), plus a `p3d` scratch pointer, and routed the
remap through it (`p3d => PLE; PLE(1:IM,1:JM,0:LM) => p3d`, etc.) -
**not** a self-remap, since gfortran-15 rejects that even on a
CONTIGUOUS pointer (confirmed via isolated test compile while fixing
the IRRAD instance of this same bug - see that file's notes). All 5 are
mandatory (no `COND`) - `PLEINST` has no conditional gating and
`FLX`/`FLC`/`FLXA`/`FLA` are forced-allocated by the
`MAPL_GridCompAddConnection(..., dst_comp="<self>")` wiring in
`SetServices` - so no `associated()` guard needed before remapping,
same reasoning as `Run`'s Edge remap block in `GEOS_IrradGridComp.F90`.

## Build verification caveat
Editing `CMakeLists.txt` forces a full top-level CMake reconfigure (not
incremental) on the next `make`. In this session's sandboxed tool
environment that reconfigure failed at `find_package(MPI REQUIRED)` inside
`src/Shared/@MAPL/superstructure/generic/CMakeLists.txt` even though the
underlying MPI toolchain (spack openmpi) was independently verified to work
fine when invoked directly - looked like a `DYLD_LIBRARY_PATH`/try-run
environment quirk specific to the non-interactive shell, not a defect in
the ported code. **The port has not been build-verified end-to-end** -
build it in a normal interactive terminal before trusting it compiles.
