# GEOSradiation_GridComp: previous architectural decisions

**Audience:** coding agents and developers working in the `GEOSradiation_GridComp` repository.
**Describes:** the MAPL3 line, branch `release/MAPL-v3`.
**Verified against:** `GEOSradiation_GridComp` @ `13f8eb8` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `f50622d`).
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

Only what is easy to see in the code now. Entries grow over time. Context:
[`ARCHITECTURE.md`](ARCHITECTURE.md).

### Full transfer at a refresh, cheap update every step
- Pattern: `SOLAR` and `IRRAD` run a full radiative transfer only when their alarm rings (`SOLAR_DT`, `IRRAD_DT`) and store normalized fluxes (solar) or fluxes plus surface-temperature derivative (infrared); every call recomputes the exports from that storage with the current insolation or surface temperature (`GEOSsolar_GridComp/GEOS_SolarGridComp.F90`, `Run`, `UPDATE_EXPORT`; `GEOSirrad_GridComp/GEOS_IrradGridComp.F90`, `Run`, `Update_Flx`). [Verified]
- Likely reason: cost of the transfer; insolation and surface temperature change fastest. [Inferred]
- Question for maintainers: heating rates from `RADIATION` do not follow the atmosphere between refreshes (`ARCHITECTURE.md`, section 5). Is that the accuracy intended at the default intervals, for both schemes RRTMG and RRTMGP?

### Independent solar and infrared scheme keys, Chou-Suarez as the default
- Pattern: four keys (`USE_RRTMGP_SORAD`, `USE_RRTMG_SORAD`, `USE_RRTMGP_IRRAD`, `USE_RRTMG_IRRAD`) are read as real numbers by shared routines; the scheme with no key set is Chou-Suarez (`GEOS_RadiationShared/rad_utils.F90`, `choose_solar_scheme`, `choose_irrad_scheme`). [Verified]
- Likely reason: the Chou-Suarez code is the original scheme, and mixed pairs were meant to be possible. [Inferred]
- Question for maintainers: are mixed pairs (for example RRTMGP shortwave with Chou-Suarez longwave) supported for production use, or only for testing?

### Band-count consistency test instead of derived `NUM_BANDS`
- Pattern: both children compute the expected total of solar plus infrared bands from the schemes and compare with the resource `NUM_BANDS`, failing the run on a mismatch (`GEOS_IrradGridComp.F90`, `Run`; `GEOS_SolarGridComp.F90`, `Run`). [Verified]
- Likely reason: the aerosol-optics tables of the provider are chosen by the same band count, so a mismatch would silently use wrong tables; a check is cheaper than a derivation across components. [Inferred]
- Question for maintainers: is `NUM_BANDS` meant to remain a hand-set key, or to be derived from the scheme choice?

### Provider-supplied aerosol optics through a state item
- Pattern: `SOLAR` and `IRRAD` import the Chemistry `AERO` state item and call the optics method of the provider for each band, instead of reading optical tables themselves (`GEOS_SolarGridComp.F90`, `Run`; `GEOS_IrradGridComp.F90`, `compute_provider_aerosol_optics` in both; the MAPL2 Physics code, not compiled on this line, makes the connection: `GEOS_PhysicsGridComp.F90`, `SetServices`, `AERO`). The header of `SOLAR` still describes an older global-table method (solar document, Observation 2). [Verified]
- Likely reason: the aerosol model owns the optical tables and the humidity growth; radiation needs only optical depth, albedo and asymmetry per band. [Inferred]
- Question for maintainers: see `ARCHITECTURE.md`, Observation 8, on the label of the optics method in MAM, DNA and the GEOS-Chem coupler.

### Shared nested RRTMGP built through a wrapper CMake
- Pattern: the `RRTMGP` repository is checked out at `GEOS_RadiationShared/@RRTMGP`, but its sources are listed and built by this repository's `GEOS_RadiationShared/RRTMGP_cmake/CMakeLists.txt`, with GEOS-specific random-number code (`rng/`) added (`components.yaml`, `RRTMGP`; `GEOS_RadiationShared/RRTMGP_cmake/CMakeLists.txt`). [Verified]
- Likely reason: the upstream build files do not match the GEOS build system or the need for MKL or Mersenne Twister generators. [Inferred]

### The composite promotes only a subset of child exports
- Pattern: `RADIATION` declares its own exports (heating rates, `ALW`, `BLW`, `RADSRF`) and re-exports a small set from `SOLAR` and `IRRAD` with `MAPL_GridCompReexport`; History reads all other child exports under the child names (`GEOS_RadiationGridComp.F90`, `SetServices`; `regression/radiation/history.yaml`; `HISTORY.AGCM.rc.tmpl` for the MAPL2 scripts). On this line the composite gets the child fields it combines through connections to `"<self>"` into its own imports, instead of requesting child exports with `alloc=.TRUE.` (see the entry below). [Verified]
- Likely reason: only the promoted fields are used by Physics siblings. [Inferred]

### Satellite simulator as an optional third child
- Pattern: `SATSIM` is added only if the sum of the `USE_SATSIM*` keys is positive, and receives `MCOSZ` from `SOLAR` by a connection inside the composite (`GEOS_RadiationGridComp.F90`, `SetServices`). [Verified]
- Likely reason: the simulators are costly diagnostics, and need the radiation grid, clouds and sun geometry. [Inferred]
- Question for maintainers: The `gcm_setup` script of the fixture (MAPL2-style) sets `USE_SATSIM_ISCCP: 1` by default, so the child is present in new experiments, and the regression cases that include `SATSIM` also set it to 1. Is that the intended default?

### Internal state items that are also exports
- Pattern: `SOLAR` marks refresh-time results with `ADD2EXPORT` in `Solar_StateSpecs.rc` (the MAPL2 line used `FRIENDLYTO = <component name>`), so that they are exports as well and are reloaded at replay; `IRRAD` uses `ADD2EXPORT` for `DFDTS` and the RATS `DFDTS_*` fields (`GEOS_SolarGridComp.F90`, `SetServices`, comment; `Solar_StateSpecs.rc`; `Irrad_StateSpecs.rc`, `DFDTS`). [Verified]
- Likely reason: the source comment in `SOLAR` (kept from the MAPL2 code) says MAPL reloads imports and internals, not exports, at the replay corrector phase; whether MAPL3 behaves so was not checked. [Verified from comment]

### Cloud-sampling settings held in shared module state
- Pattern: the composite `Initialize` sets the condensate inhomogeneity and the correlation lengths in module variables of `GEOS_RadiationShared` (`set_inhomogeneity`, `initialize_cloud_subcol_gen`), which both children read when they sample subcolumns (`GEOS_RadiationGridComp.F90`, `Initialize`). On this line `IRRAD` and `SOLAR` also call `set_inhomogeneity` in their own `Initialize`, with a comment that this is for standalone runs and to be removed (`GEOS_IrradGridComp.F90`, `GEOS_SolarGridComp.F90`, `Initialize`). [Verified]
- Likely reason: one setting for both bands, read once, without passing through the import and export states. [Inferred]
- Question for maintainers: do the two children share these settings by intent? With the module state they cannot differ at present.

### RRTMGP wrapper state held on the component
- Pattern: each child attaches a private state (`ty_RRTMGP_state`) to the gridded component (on this line with `_SET_NAMED_PRIVATE_STATE` in `SetServices`) and loads the k-distribution tables on first use (`GEOS_SolarGridComp.F90`, `SetServices`; `GEOS_IrradGridComp.F90`, `SetServices`, `LW_Driver`). [Verified]
- Likely reason: avoid reading the tables at every refresh and avoid a global variable. [Inferred]

### Composite reads child fields through connections to itself (MAPL3 line)
- Pattern: `RADIATION` declares the child fluxes it needs as its own imports (`Radiation_StateSpecs.rc`) and fills them with `MAPL_GridCompAddConnection(..., dst_comp="<self>")`; `Run` runs the children one at a time with `MAPL_GridCompRunChild` and then reads its imports (`GEOS_RadiationGridComp.F90`, `SetServices`, `Run`; comments there). [Verified]
- Likely reason: MAPL3 has no direct way to read a child's export state in the parent; the porting note says so (`mapl3-porting-notes.md`, "Parent Run reading a child's export state"). [Verified from note]

### Each child creates its own refresh alarm (MAPL3 line)
- Pattern: `IRRAD`, `SOLAR` and `SATSIM` create a named alarm in `Initialize` from `<name>_DT`, with the code of MAPL2's `handle_clock_and_main_alarm` copied in (`GEOS_IrradGridComp.F90`, `GEOS_SolarGridComp.F90`, `GEOS_SatsimGridComp.F90`, `Initialize`; comments there). [Verified]
- Likely reason: MAPL3 no longer creates the alarm for a component (stated in those comments). [Verified from comment]
- Question for maintainers: three copies of the same alarm set-up exist; is a shared routine wanted?

### SatSim masked exports removed in the port (MAPL3 line)
- Pattern: the optional `Masked_Exports::` table, the `SATORB` import and the masking pass were dropped; the code comment says the removed code was last available, commented out, at commit `0a0936c`, and that the live MAPL2 version is on `release/MAPL-v3` at `7862c74` (`GEOS_SatsimGridComp.F90`, `SetServices`, comment; the commit `7862c74` is listed on this branch as a Solar change, so the comment's pointer does not match it). [Verified for the code comment; *Unverified* for the `7862c74` pointer]
- Likely reason: the comment says they were MAPL2-only and never used. [Verified from comment]
