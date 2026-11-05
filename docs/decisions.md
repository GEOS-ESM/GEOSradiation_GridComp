# GEOSradiation_GridComp: previous architectural decisions

**Audience:** coding agents and developers working in the `GEOSradiation_GridComp` repository.
**Describes:** the MAPL2 line, branch `develop`.
**Verified against:** `GEOSradiation_GridComp` @ `0b402f6` (2026-11-05).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `0f2e0fb`).
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

Only what is easy to see in the code now. Entries grow over time. Context:
[`ARCHITECTURE.md`](ARCHITECTURE.md).

### Full transfer at a refresh, cheap update every step
- Pattern: `SOLAR` and `IRRAD` run a full radiative transfer only when their alarm rings (`SOLAR_DT:`, `IRRAD_DT:`) and store normalized fluxes (solar) or fluxes plus surface-temperature derivative (infrared); every call recomputes the exports from that storage with the current insolation or surface temperature (`GEOSsolar_GridComp/GEOS_SolarGridComp.F90`, `Run`, `UPDATE_EXPORT`; `GEOSirrad_GridComp/GEOS_IrradGridComp.F90`, `Run`, `Update_Flx`). [Verified]
- Likely reason: cost of the transfer; insolation and surface temperature change fastest. [Inferred]
- Question for maintainers: heating rates from `RADIATION` do not follow the atmosphere between refreshes (parent `ARCHITECTURE.md`, section 5). Is that the accuracy intended at the default intervals, for both schemes RRTMG and RRTMGP?

### Independent solar and infrared scheme keys, Chou-Suarez as the default
- Pattern: four keys (`USE_RRTMGP_SORAD:`, `USE_RRTMG_SORAD:`, `USE_RRTMGP_IRRAD:`, `USE_RRTMG_IRRAD:`) are read as real numbers by shared routines; the scheme with no key set is Chou-Suarez (`GEOS_RadiationShared/rad_utils.F90`, `choose_solar_scheme`, `choose_irrad_scheme`). [Verified]
- Likely reason: the Chou-Suarez code is the original scheme, and mixed pairs were meant to be possible. [Inferred]
- Question for maintainers: are mixed pairs (for example RRTMGP shortwave with Chou-Suarez longwave) supported for production use, or only for testing?

### Band-count consistency test instead of derived `NUM_BANDS`
- Pattern: both children compute the expected total of solar plus infrared bands from the schemes and compare with the resource `NUM_BANDS:`, failing the run on a mismatch (`GEOS_IrradGridComp.F90`, `Run`; `GEOS_SolarGridComp.F90`, `Run`). [Verified]
- Likely reason: the aerosol-optics tables of the provider are chosen by the same band count, so a mismatch would silently use wrong tables; a check is cheaper than a derivation across components. [Inferred]
- Question for maintainers: is `NUM_BANDS:` meant to remain a hand-set key, or to be derived from the scheme choice?

### Provider-supplied aerosol optics through a state item
- Pattern: `SOLAR` and `IRRAD` import the Chemistry `AERO` state item and call the optics method of the provider for each band, instead of reading optical tables themselves (`GEOS_SolarGridComp.F90`, `Run`; `GEOS_IrradGridComp.F90`, `compute_provider_aerosol_optics`; `GEOS_PhysicsGridComp.F90`, `SetServices`, `AERO`). The header of `SOLAR` and the `Run` description of the composite still describe an older global-table method (parent `ARCHITECTURE.md`, Observation 1; solar document, Observation 2). [Verified]
- Likely reason: the aerosol model owns the optical tables and the humidity growth; radiation needs only optical depth, albedo and asymmetry per band. [Inferred]
- Question for maintainers: see parent `ARCHITECTURE.md`, Observation 8, on the label of the optics method in MAM, DNA and the GEOS-Chem coupler.

### Shared nested RRTMGP built through a wrapper CMake
- Pattern: the `RRTMGP` repository is checked out at `GEOS_RadiationShared/@RRTMGP`, but its sources are listed and built by this repository's `GEOS_RadiationShared/RRTMGP_cmake/CMakeLists.txt`, with GEOS-specific random-number code (`rng/`) added (`components.yaml`, `RRTMGP`; `GEOS_RadiationShared/RRTMGP_cmake/CMakeLists.txt`). [Verified]
- Likely reason: the upstream build files do not match the GEOS build system or the need for MKL or Mersenne Twister generators. [Inferred]

### The composite promotes only a subset of child exports
- Pattern: `RADIATION` declares its own exports (heating rates, `ALW`, `BLW`, `RADSRF`) and promotes a small set from `SOLAR` and `IRRAD`; History reads all other child exports under the child names (`GEOS_RadiationGridComp.F90`, `SetServices`, `Run`; `HISTORY.AGCM.rc.tmpl`). The composite requests the child exports it needs with `alloc=.TRUE.` before running the children. [Verified]
- Likely reason: only the promoted fields are used by Physics siblings. [Inferred]

### Satellite simulator as an optional third child
- Pattern: `SATSIM` is added only if the sum of the `USE_SATSIM*` keys is positive, and receives `MCOSZ` from `SOLAR` by a connection inside the composite (`GEOS_RadiationGridComp.F90`, `SetServices`). [Verified]
- Likely reason: the simulators are costly diagnostics, and need the radiation grid, clouds and sun geometry. [Inferred]
- Question for maintainers: `gcm_setup` sets `USE_SATSIM_ISCCP: 1` by default, so the child is present in new experiments. Is that the intended default?

### Internal state items that are also exports
- Pattern: `SOLAR` marks refresh-time results with `FRIENDLYTO = <component name>`, so that they are exports as well and are reloaded at replay; `IRRAD` uses `ADD2EXPORT` for `DFDTS` and the RATS `DFDTS_*` fields (`GEOS_SolarGridComp.F90`, `SetServices`; `Irrad_StateSpecs.rc`, `DFDTS`). [Verified]
- Likely reason: the source comment in `SOLAR` says MAPL reloads imports and internals, not exports, at the replay corrector phase. [Verified from comment]

### Cloud-sampling settings held in shared module state
- Pattern: the composite `Initialize` sets the condensate inhomogeneity and the correlation lengths in module variables of `GEOS_RadiationShared` (`set_inhomogeneity`, `initialize_cloud_subcol_gen`), which both children read when they sample subcolumns (`GEOS_RadiationGridComp.F90`, `Initialize`). [Verified]
- Likely reason: one setting for both bands, read once, without passing through the import and export states. [Inferred]
- Question for maintainers: do the two children share these settings by intent? With the module state they cannot differ at present.

### RRTMGP wrapper state held on the component
- Pattern: each child attaches a private state (`ty_RRTMGP_state`) to the gridded component with `ESMF_UserCompSetInternalState` and loads the k-distribution tables on first use (`GEOS_SolarGridComp.F90`, `SetServices`; `GEOS_IrradGridComp.F90`, `SetServices`, `LW_Driver`). [Verified]
- Likely reason: avoid reading the tables at every refresh and avoid a global variable. [Inferred]
