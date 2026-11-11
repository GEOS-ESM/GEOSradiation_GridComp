# GEOSradiation_GridComp architecture

A map of the repository: what lives where, how the pieces fit, and what must not break. Paths are relative to the repository root. It names modules and routines rather than linking to lines; use them to search the code. It describes only things that change rarely; status, rules for changes, known traps and to-dos are in [`AGENTS.md`](../AGENTS.md).

## Bird's-eye view

Radiation computes the solar (shortwave) and infrared (longwave) radiative fluxes of the GEOS atmosphere and the heating rates and surface terms derived from them, plus optional satellite-simulator diagnostics. In GEOS it is the `RADIATION` child of Physics, run last, after the components whose clouds, gases, aerosols and surface it uses. It is the `RADIATION` composite with three children, `SOLAR`, `IRRAD` and `SATSIM`, and code they share.

The defining design choice: a full radiative transfer is expensive, so each child does it only at its refresh interval and, on every other step, updates its fluxes cheaply from what the last refresh stored (see "Refresh and every-step update").

## Codemap

- `GEOS_RadiationGridComp.F90`, `Radiation_StateSpecs.rc`: the `RADIATION` composite. Adds the children, runs them, and combines their fluxes into heating rates and surface terms. No internal state.
- `GEOSsolar_GridComp/` (`GEOS_SolarGridComp.F90`, `Solar_StateSpecs.rc`): `SOLAR`, shortwave. Schemes Chou-Suarez (`sorad.F90`), RRTMG (`RRTMG/rrtmg_sw`) and RRTMGP; load balancing between day and night columns (`SORADCORE`).
- `GEOSirrad_GridComp/` (`GEOS_IrradGridComp.F90`, `Irrad_StateSpecs.rc`): `IRRAD`, longwave. Schemes Chou-Suarez (`irrad.F90`), RRTMG (`RRTMG/rrtmg_lw`) and RRTMGP; optional per-tracer diagnostics (`RATS_DIAGNOSTICS`).
- `GEOSsatsim_GridComp/` (`GEOS_SatsimGridComp.F90`, `Satsim_StateSpecs.rc`): `SATSIM`, diagnostic satellite simulators from COSP (ISCCP, MODIS, MISR, radar in `quickbeam/`, lidar in `actsim/`). Feeds no tendency.
- `GEOS_RadiationShared/`: code both bands use: scheme choice (`rad_utils.F90`), stochastic cloud subcolumns (`cloud_subcol_gen.F90`), condensate inhomogeneity (`cloud_condensate_inhomogeneity.F90`), Chou-Suarez cloud optics (`gettau.F90`, `radconstants.F90`); the CMake wrapper and random-number code for RRTMGP (`RRTMGP_cmake/`); scripts that prepare the NRLSSI2 solar-variability file (`NRLSSI2/`).
- `GEOS_RadiationShared/@RRTMGP`: RTE-RRTMGP, a separate nested repository, compiled through `RRTMGP_cmake/`.
- `regression/`: standalone regression cases, and `FakeGOCART_GridComp`, an aerosol stand-in that provides an empty `AERO` state.
- `mapl3-porting-notes.md` (root and in each child directory): history of the MAPL3 ports.

## How Radiation is put together

- **Children.** `SetServices` adds `IRRAD` and `SOLAR` always and `SATSIM` when any `USE_SATSIM*` key is positive (`get_use_satsim`), each with its own configuration file (`irrad.yaml`, `solar.yaml`, `satsim.yaml`).
- **Run order.** The composite's single `run` phase runs `IRRAD`, `SOLAR`, then `SATSIM` (which needs `SOLAR`'s `MCOSZ`), then combines: `DTDT` and the per-band heating rates (`RADLW`, `RADSW`, their clear-sky and no-aerosol variants) from flux divergence, `RADSRF`, and the linearized surface longwave `ALW`, `BLW`.
- **The composite reads its children's exports through itself.** Fields it combines are declared as its own imports and filled by connections from the children to `"<self>"` (`FLX`, `FLC`, `FLXA`, `FLA`, `DSFDTS0`, `SFCEM0`, `TSREFF` from `IRRAD`; `FSW`, `FSC`, `FSWNA`, `FSCNA` from `SOLAR`). A few child exports are re-exported as the composite's own (from `SOLAR`: the PAR, NIR and UV surface fluxes, cloud optical depths and fractions, `ALBEDO`, band fluxes; from `IRRAD`: `LWS`, `LWS0`).
- **Shared settings.** The composite's `Initialize` sets module-level settings that both bands use: condensate inhomogeneity (`RAD_CONDENSATE_INHOMOGENEITY`) and the cloud-overlap correlation lengths (`ADL_*`, `RDL_*`).

### Refresh and every-step update

Each child creates its own alarm in `Initialize` (`solar_run_alarm` from `SOLAR_DT`, `irrad_lw_alarm` from `IRRAD_DT`, `satsim_alarm` from `SATSIM_DT`; reference time from `REFERENCE_DATE`, `REFERENCE_TIME`, `RUN_AT_INTERVAL_START`).

| | Refresh (alarm rings) | Every step |
|---|---|---|
| `SOLAR` | full transfer for the mean sun position over the next interval; fluxes normalized by top-of-atmosphere insolation, stored as internal state | stored normalized fluxes times current insolation (`UPDATE_EXPORT`) |
| `IRRAD` | full transfer; fluxes and their derivative with respect to surface temperature, stored as internal state | stored fluxes plus derivative times the change in surface temperature (`Update_Flx`) |
| `SATSIM` | everything | nothing |

`CALLED_LAST` sets whether the update comes before the refresh (1, the default) or after it (0).

### Scheme selection

Solar and infrared schemes are chosen independently by `USE_RRTMGP_SORAD`, `USE_RRTMG_SORAD`, `USE_RRTMGP_IRRAD`, `USE_RRTMG_IRRAD` (non-zero means on), with RRTMGP before RRTMG before Chou-Suarez, the default (`choose_solar_scheme`, `choose_irrad_scheme` in `rad_utils.F90`). Band counts: solar 8 (Chou-Suarez) or 14 (RRTMG, RRTMGP); infrared 10 or 16. `NUM_BANDS` must equal the solar plus infrared bands of the chosen pair.

## Contracts with the rest of GEOS

### Imports

- From Moist: water vapour and condensates (`QV`, `QL`, `QI`, `QR`, `QS`, `QG`), cloud fraction `FCLD`, effective radii (`RL`, `RI`, `RR`, `RS`, `RG`).
- From Surface: albedos (`ALBVR`, `ALBVF`, `ALBNR`, `ALBNF`), emissivity `EMIS`, skin temperature `TS` (also as `TSINST`), land and ocean fractions.
- From Chemistry: radiatively active gases (`O3`, `CH4`, `N2O`, `CFC11`, `CFC12`, `HCFC22`, optionally a 3-D `CO2`) and the `AERO` state for aerosol optics (Appendix).
- From the Agcm and Physics: `PLE`, `PLEINST`, `T`, `PREF`, `ZLE`.

### Exports

- `DTDT` (pressure-weighted heating) to the Physics tendencies; `RADLW` and `RADSW` to Moist; `RADLW` and its clear-sky version `RADLWC` to Turbulence.
- To Surface: `ALW`, `BLW`, the normalized PAR, NIR and UV surface fluxes, `LWS0` (as downward longwave at the surface), and `DROBIO`, `DFOBIO` when `USE_OCEANOBIOGEOCHEM` is on.
- To Chemistry: cloud optical depths and fraction (`TAUCLI`, `TAUCLW`, `CLDTT`), PAR (`DFPAR`, `DRPAR`).
- Every other child export (band fluxes, diagnostics, all of `SATSIM`'s) is available only under the child's name.

### Services

Radiation neither provides nor subscribes to any service.

## Invariants

- Fluxes are net downward on interfaces indexed `0:LM` from the top of the atmosphere; heating is top minus bottom. MAPL3 gives Edge fields bounds `1:LM+1`, so every routine that uses the `0:LM` convention re-bounds its Edge pointers first.
- `NUM_BANDS` equals the solar plus infrared band counts of the chosen schemes, and bands form one sequence: solar first, infrared after (`BANDS_SOLAR_OFFSET`, `OFFSET`). The aerosol provider's band numbering follows it.
- Between refreshes, fluxes change only through insolation (shortwave) and surface temperature (longwave); new clouds, humidity or albedo reach the heating at the next refresh.
- Refresh results are internal state; a restart continues from the stored refresh, so restarts are seamless only at refresh times.
- Condensate inhomogeneity and correlation lengths are module-level and shared by both bands, not per child.

## Cross-cutting concerns

- **Generated specs:** each component's imports, exports and internals are declared in its `*_StateSpecs.rc` and turned into `*___.h` headers by MAPL's ACG. Added by hand in `SetServices`, because the table cannot express them: the `AERO` state import, `IRRAD`'s `RATS_DIAGNOSTICS` fields, and `SOLAR`'s per-band internals and exports, built in loops.
- **Exports only when requested:** most exports are computed only if requested, and some work is too (`SOLAR`'s no-aerosol transfer runs only when a no-aerosol export is requested; `IRRAD`'s choice of flux sets under RRTMGP).
- **Undefined values:** `SOLAR`, `IRRAD` and `SATSIM` use MAPL2's undefined value (`MAPL2_UNDEF = 1.0e15`), not MAPL3's `MAPL_UNDEF`; for example `SATSIM`'s exports and RRTMG's no-aerosol infrared fluxes.
- **Build options:** `ENABLE_SOLAR_RADVAL` adds developer diagnostics to `SOLAR` and changes its restart contents; `RRTMGP_SINGLE_PRECISION` builds RRTMGP in single precision; RRTMGP's random numbers come from MKL when it is found, otherwise from a Mersenne Twister.
- **Testing:** each regression case (`regression/`) runs `IRRAD`, `SOLAR`, `SATSIM` or the whole composite for one step under a configurable Root, with every import supplied by a data component restarted from a MAPL2 capture and an empty `AERO` from `FakeGOCART_GridComp`, and compares the checkpoints exactly against a stored MAPL3 baseline that was validated against the MAPL2 capture.

## Where to go next

- [`AGENTS.md`](../AGENTS.md): rules, status on this line, to-dos, known traps, build and test.
- `mapl3-porting-notes.md` at the root and in `GEOSsolar_GridComp/`, `GEOSirrad_GridComp/`, `GEOSsatsim_GridComp/`: how the ports were done.

## Appendix: aerosol optics from the `AERO` state

`SOLAR` and `IRRAD` get aerosol optics from the provider's `AERO` state (in GEOS, the Chemistry component named by `AERO_PROVIDER`, normally GOCART2G) at each refresh, through a callback the provider attaches to that state (`compute_provider_aerosol_optics` in each child):

1. If the state's Info key `implements_aerosol_optics_method` is false, no aerosols are used.
2. Write relative humidity and interface pressure into the provider fields named by `relative_humidity_for_aerosol_optics` and `air_pressure_for_aerosol_optics`.
3. For each band, set `band_for_aerosol_optics`, run the callback `run_aerosol_optics` (`ESMF_MethodExecute`, checking both `rc` and `userRC`), and read the fields named by `extinction_in_air_due_to_ambient_aerosol`, `single_scattering_albedo_of_ambient_aerosol` and `asymmetry_parameter_of_ambient_aerosol`. They hold extinction, extinction times albedo, and extinction times albedo times asymmetry; `SOLAR` divides to get albedo and asymmetry. Both children clip these fields as received, before any division: all three at zero below, the second at 0.999999 and the third at 0.999 above (as in MAPL2).
