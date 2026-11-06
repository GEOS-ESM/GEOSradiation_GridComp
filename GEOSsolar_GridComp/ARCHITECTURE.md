# GEOSsolar_GridComp (SOLAR): architecture

**Audience:** coding agents and developers working on the shortwave (solar) radiation gridded component.
**Describes:** the MAPL3 line, branch `release/MAPL-v3`.
**Verified against:** `GEOSradiation_GridComp` @ `13f8eb8` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `f50622d`; repository `GEOSradiation_GridComp`, mounted at `src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSradiation_GridComp`): `.../@GEOSradiation_GridComp/GEOSsolar_GridComp`.
Other fixtures (for example GEOSadas, GEOSldas) may place this repository elsewhere; paths inside this repository do not change.

Paths are relative to `GEOSsolar_GridComp/` unless they start with `GEOS_RadiationShared/` (same repository root) or `src/` (fixture-relative).
Parent: [Radiation](../docs/ARCHITECTURE.md) (scheme choice, refresh versus update, aerosol optics protocol, interfaces to Physics siblings: read it first). Sibling: [Irrad](../GEOSirrad_GridComp/ARCHITECTURE.md). Repository index: [../AGENTS.md](../AGENTS.md). Porting history: `mapl3-porting-notes.md` in this directory (not repeated here). Build status on this line: built as part of the `RADIATION` library and run in the `solar-sa` and `radiation` regression cases (`../CMakeLists.txt`, `../regression/cases.txt`).

## 1. Purpose
`SOLAR` computes shortwave (solar) radiative fluxes in a cloudy, aerosol-laden atmosphere with one of three schemes: Chou-Suarez `SORAD` (default), RRTMG, or RRTMGP. A full transfer calculation runs only when the component's alarm rings; every call fills the exports from stored normalized fluxes and the current sun position (`GEOS_SolarGridComp.F90`, header comment, `Run`).

## 2. Place in the hierarchy
- Parent: `RADIATION`, which adds it as `SOLAR` with `MAPL_GridCompAddChild(gc, "SOLAR", solarSetServices, "solar.yaml", ...)` (`GEOS_RadiationGridComp.F90`, `SetServices`). It has no children. The file ends with a routine `Solar_SetServices` that wraps `SetServices` as a library entry point (the `solar_setservices_` of `../regression/solar-sa/root.yaml`).
- Files: `GEOS_SolarGridComp.F90` (module `GEOS_SolarGridCompMod`, about 6890 lines, everything below); `Solar_StateSpecs.rc` (ACG input); `sorad.F90` (module `soradmod`, Chou-Suarez `SORAD`); `soradconstants.F90` (module `sorad_constants`); `RRTMG/rrtmg_sw/gcm_model/` (AER RRTMG shortwave: `modules/` data, `src/` code, `NRLSSI2.F90` solar-variability model; libraries `RRTMG_SW`, `RRTMG_SW_mods`); `SolarDoc.pdf` (Chou-Suarez report, referenced in the module header).
- RRTMGP code is not in this directory: it is the nested `RRTMGP` repository built by `GEOS_RadiationShared/RRTMGP_cmake` (parent document, section 1).

## 3. Entry points
- `SetServices` (`GEOS_SolarGridComp.F90`, `SetServices`): calls `choose_solar_scheme` to set the number of solar bands (the ungridded dimension `num_bands_solar`); reads `USE_OCEANOBIOGEOCHEM` (adds the ungridded dimension `nb_obio`, 33 bands, and sets the condition `SOLAR_TO_OBIO`); sets `USE_SOLAR_RADVAL` from the compile option; includes the ACG-generated import, internal and export specs from `Solar_StateSpecs.rc`; adds by hand the state-item import `AERO` (`MAPL_GridCompAddSpec`, `itemtype=MAPL_STATEITEM_STATE`) and the per-band `OSRBbbRGN`, `ISRBbbRGN` internals and `OSRBbbRG`, `ISRBbbRG`, `TBRBbbRG` exports in loops (RRTMG and RRTMGP only); registers `Initialize` and `Run` (`phase_name='run'`). Finalize is the MAPL default. The routine has no `MAPL_GenericSetServices` call and does not read `RUN_DT`.
- `Initialize` (`GEOS_SolarGridComp.F90`, `Initialize`): creates the alarm `solar_run_alarm` with interval `<component name>_DT` (`SOLAR_DT`, default the time step) with the alarm set-up of MAPL2's generic code copied in (parent document, section 5); creates the solar orbit with `MAPL_SunOrbitCreate` from `ECCENTRICITY`, `OBLIQUITY`, `PERIHELION`, `EQUINOX`, `EOT`, `ORBIT_ANAL2B` and the `ORB2B_*` keys (defaults in the code) and stores it in the private state, because MAPL3 does not provide an orbit to the component; calls `set_inhomogeneity` (comment: for standalone runs, to be removed).
- `Run`, one phase (`GEOS_SolarGridComp.F90`, `Run`):
  1. `PRELIMS`: fetch the alarm and the orbit; read `PRS_LOW_MID_CLOUDS` (default 70000 Pa), `PRS_MID_HIGH_CLOUDS` (40000 Pa), `CO2`, `SOLAR_CONSTANT`, `SUN_FLAG` (0), `DYCORE`, `SOLAR_LOAD_BALANCE` (1), `SOLAR_LB_MAX_PASSES` (100), `NUM_BANDS`, `CALLED_LAST` (1). `CO2`, `SOLAR_CONSTANT`, `DYCORE` and `NUM_BANDS` have no default. Time-dependent `CO2` (-1) calls `getco2`; scheme choice and the `NUM_BANDS` consistency test (parent document, section 4); solar-constant source (below); find the model levels `LCLDMH`, `LCLDLM` that split high, middle and low clouds from the import `PREF`.
  2. `UPDATE_EXPORT` first if `CALLED_LAST /= 0`.
  3. `REFRESH` if `ESMF_AlarmIsRinging` for `solar_run_alarm`: optical properties of the aerosol provider (`compute_provider_aerosol_optics`, parent document, section 6); a without-aerosol `SORADCORE` call only if a no-aerosol export is requested, else zero the no-aerosol internals; the regular `SORADCORE` with aerosols. The refresh covers the mean sun position over the alarm ring interval starting at the current time plus one clock step when the update came first, or at the current time when it came last.
  4. `UPDATE_EXPORT` after the refresh if `CALLED_LAST == 0`.
- `SORADCORE` (contained in `Run`), per refresh: gets the interval-mean insolation from `MAPL_SunGetInsolation` (`INTV` = alarm ring interval); packs the sunlit columns of all imports into a buffer; optionally redistributes the work among processors (`MAPL_LoadBalancerCreate`, `MAPL_LoadBalancerRun`, `MAPL_LoadBalancerDestroy`; the load-balancer interface was renamed in MAPL, see the recent commits of this branch); calls the chosen scheme; retrieves; unpacks into the internal state, filling night points with the spec default; stores fluxes divided by the absolute insolation (normalized). Scheme branch `SCHEME`: `USE_CHOU` calls `shrtwave` (contained, wraps `SORAD`); `USE_RRTMGP` uses `PROCESS_RRTMGP_BLOCK` per block of columns in an OpenMP loop; `USE_RRTMG` calls `RRTMG_SW`.
- `UPDATE_EXPORT` (contained in `Run`), every call: `MAPL_SunGetInsolation` for the current clock step (`ZTH`, `SLR`, scaled by the solar constant), then each flux export = normalized internal flux times `SLR` (surface band fluxes, `FSW`, `FSC`, `FSWU`, `FSWD`, per-band `FSWBAND`, `OSR`, `RSR`, `SLRSF` and their clear-sky and no-aerosol versions); surface albedo diagnostic `ALBEDO` from the normalized surface fluxes; `MCOSZ` = mean cosine of the zenith angle over the step, `COSZ`; cloud diagnostics.

### Cloud diagnostics at the heartbeat
When any of `TAUCLI`, `TAUCLW`, `TAUCLR`, `TAUCLS`, `TAULO`, `TAUMD`, `TAUHI`, `COTxx`, `CLDTMP`, `CLDPRS` is requested, `UPDATE_EXPORT` calls `GETVISTAU` (module `gettau`) column by column on the current imports. Chemistry consumes `TAUCLI`, `TAUCLW` (parent document, section 7), so these are refreshed every step, not every refresh. `CLDHI`, `CLDMD`, `CLDLO`, `CLDTT` are derived from the imported `FCLD`. The `...SWHB` fields (subcolumn cloud fractions at the heartbeat) exist only with `ENABLE_SOLAR_RADVAL`.

## 4. State
- Imports (`Solar_StateSpecs.rc`, category IMPORT): `PLE`, `TS`, `CH4`, `N2O`, `T`, `QV`, `QL`, `QI`, `QR`, `QS`, `QG`, `RL`, `RI`, `RR`, `RS`, `RG`, `O3` (mass mixing ratio, `kg kg-1`), `FCLD`, surface albedos `ALBVR`, `ALBVF`, `ALBNR`, `ALBNF` and reference pressure `PREF`; state item `AERO` added by hand (`itemtype = MAPL_STATEITEM_STATE`). The specs have no averaging or refresh interval on this line. Who provides them: parent document, section 7 (in the regression cases, a data component: `../regression/solar-sa/root.yaml`). `CH4`, `N2O` and the surface `TS` are declared but this component's shortwave gas set uses `h2o`, `co2`, `o3`, `n2o`, `co`, `ch4`, `o2`, `n2` for RRTMGP (`SORADCORE`, `gas_concs%init`).
- Internal state (`Solar_StateSpecs.rc`, category INTERNAL, about 270 rows in the file in all categories, plus the per-band loops in `SetServices`): normalized net and upward fluxes `FSWN`, `FSCN`, `FSWUN`, `FSCUN` and the no-aerosol versions `FSWNAN`, `FSCNAN`, `FSWUNAN`, `FSCUNAN`; normalized surface fluxes `DRUVRN`, `DFUVRN`, `DRPARN`, `DFPARN`, `DRNIRN`, `DFNIRN`; per-band `FSWBANDN`, `FSWBANDNAN`, and `DRBANDN`, `DFBANDN` when `USE_OCEANOBIOGEOCHEM` /= 0; per-band TOA `OSRBbbRGN`, `ISRBbbRGN` for bands in `band_output_supported` (RRTMG and RRTMGP); refresh-time diagnostics `COSZSW`, `CLDxxSW`, `TAUxxPAR`, `COTxxPAR`... (the `TAUxxPAR` family only under `SOLAR_RADVAL`).
- Internals with `ADD2EXPORT` in `Solar_StateSpecs.rc` are automatically also exports of the same name, with no copy (`SetServices`, comment "INTERNALs can be EXPORTed without an explicit EXPORT"; the comment says the MAPL2 `FRIENDLYTO = trim(COMP_NAME)` is expressed by the `ADD2EXPORT` column on this line). Reason given in the source: MAPL reloads imports and internals, not exports, at the replay corrector phase, so refresh-time exports must be internals to be reloaded in step (`SetServices`; `Solar_StateSpecs.rc`, comment "The following ADD2EXPORT INTERNALs ... are really EXPORTs"). Night points take the `FILL` value of the internal row.
- Exports (`Solar_StateSpecs.rc`, category EXPORT, plus the band loops in `SetServices`); the per-band shapes use the ungridded dimensions `num_bands_solar` (14 or 8) and `nb_obio` (33), with the condition `SOLAR_TO_OBIO` for the ocean-biogeochemistry exports. The Radiation composite re-exports only a subset (parent document, section 3.3).
- Restarts: in the MAPL2-style scripts the internal state is written to `solar_internal_rst` (`@MODELATMSOLAR_INTERNAL_RESTART_FILE` in `src/Applications/@GEOSgcm_App/AGCM.rc.tmpl`); the regression cases restart and checkpoint it (`../regression/solar-sa/solar.yaml`). `AERO` is not an internal. The RRTMGP lookup tables are not restarted; they are loaded from the files on the first refresh (below).
- RRTMGP private state: type `ty_RRTMGP_state` (attached with `_SET_NAMED_PRIVATE_STATE`) holds `k_dist`, an `initialized` flag and the solar orbit; it is loaded once from `RRTMGP_GAS_SW` (default `rrtmgp-gas-sw-g112.nc`) on the first call that needs it, which can be `SORADCORE` or `UPDATE_EXPORT` (band wavenumbers) (`SORADCORE`, `UPDATE_EXPORT`). Cloud optics are loaded from `RRTMGP_CLOUD_OPTICS_SW` (`rrtmgp-clouds-sw.nc`) according to `RRTMGP_CLOUD_OPTICS_TYPE_SW` (`LUT` default, or `PADE`) on every refresh (`SORADCORE`).

## 5. Configuration that changes structure
| Key | Effect |
|---|---|
| `USE_RRTMGP_SORAD`, `USE_RRTMG_SORAD` | scheme; none set means Chou-Suarez (`GEOS_RadiationShared/rad_utils.F90`, `choose_solar_scheme`) |
| `NUM_BANDS` | must equal solar plus infrared band count |
| `SOLAR_DT`, `CALLED_LAST` | refresh interval (alarm created in `Initialize`; default the time step) and call order (update first if non-zero). `SOLARAvrg` is not read on this line |
| `USE_OCEANOBIOGEOCHEM` | adds `DROBIO`, `DFOBIO`, `DRBANDN`, `DFBANDN` |
| `SOLAR_LOAD_BALANCE`, `SOLAR_LB_MAX_PASSES`; `DYCORE: DATMO` | work redistribution among processors; `DATMO` (single column) turns it off and forces a minimum `ZTH` so that the balancer sees lit points |
| `SOLAR_CONSTANT` | `-1`: time-dependent value from `MAPL_SunGetInsolation`/`MAPL_SunGetSolarConstant` |
| `SOLAR_CYCLE_FILE_NAME` (default `/dev/null`), `USE_NRLSSI2` (true), `PERSIST_SOLAR` (true) | with a file and `USE_NRLSSI2`: NRLSSI2 forcing (solar constant, Mg and sunspot indices); allowed only for RRTMG and RRTMGP; with a file and `USE_NRLSSI2` false: band fractions for Chou-Suarez |
| `ISOLVAR` (0), `INDSOLVAR_1`, `INDSOLVAR_2`, `SOLCYCFRAC` | RRTMG solar-variability options; `ISOLVAR == 1` fails; `ISOLVAR == 2` uses the file's Mg and sunspot indices |
| `RRTMGP_GAS_SW`, `RRTMGP_CLOUD_OPTICS_SW`, `RRTMGP_CLOUD_OPTICS_TYPE_SW`, `RRTMGP_ICE_ROUGHNESS_SW`, `RRTMGP_CLOUD_OVERLAP_TYPE_SW` (default `GEN_MAX_RAN_OVERLAP`), `RRTMGP_DELTA_SCALE` (true), `RRTMGP_USE_RRTMG_ICEFLG3_LIKE_FORWICE` (true), `RRTMGP_SW_BLOCKSIZE` (4) | RRTMGP data, cloud treatment and blocking |
| `RRTMG_ICEFLG` (3), `RRTMG_LIQFLG` (1), `RRTMGSW_PARTITION_SIZE`, `RRTMGSW_USE_PRECIP_IN_RADIATION` (true if `LM > 72`) | RRTMG cloud treatment |
| `RADSW_BINARY_CLOUDS` | force cloud fraction to 0 or 1 |
| `PRS_LOW_MID_CLOUDS`, `PRS_MID_HIGH_CLOUDS`, `TAUCRIT` (0.10) | cloud-layer split; threshold for `CLDTMP`, `CLDPRS` |
| `ECCENTRICITY`, `OBLIQUITY`, `PERIHELION`, `EQUINOX`, `EOT`, `ORBIT_ANAL2B`, `ORB2B_*` | orbit, read in `Initialize` (defaults in the code) |
| `RUN_AT_INTERVAL_START`, `REFERENCE_DATE`, `REFERENCE_TIME` | alarm timing (`Initialize`) |
| Compile time: `ENABLE_SOLAR_RADVAL` | defines `SOLAR_RADVAL`; extra diagnostics (`GEOSsolar_GridComp/CMakeLists.txt`) |

Solar-forcing selection in `Run`: file with `USE_NRLSSI2` (RRTMG[P] only), else file without, else `SC < 0` (`MAPL_SunGetSolarConstant` by date), else fixed table values (`GEOS_SolarGridComp.F90`, `Run`, "Decide how to do solar forcing").

## 6. Dependencies
- Libraries (`CMakeLists.txt`): `MAPL`, `GEOS_Shared` (for `getco2`), `GEOS_RadiationShared`, `RRTMG_SW`, `RRTMG_SW_mods`, `RRTMGP`; subdirectory `RRTMG`; `HAVE_MKL` is defined when MKL is found.
- Modules used: `soradmod`, `sorad_constants`, `gettau` (`getvistau`), `rrtmg_sw_rad`, `rrtmg_sw_init`, `parrrsw`, `rrsw_wvn`, `rad_types`, `cloud_subcol_gen`, `cloud_condensate_inhomogeneity`, `rad_utils`, RRTMGP modules (`mo_gas_optics_rrtmgp`, `mo_rte_sw`, `mo_cloud_optics_rrtmgp`, `mo_fluxes_byband`, ...), `gFTL_StringVector`; MKL random numbers (`mo_rng_mklvsl_plus`) when MKL is found, else `mo_rng_mt19937`.
- Mixed-scheme pairs are not rejected by the code (for example RRTMGP shortwave with Chou-Suarez infrared) as long as `NUM_BANDS` matches; the RRTMGP and RRTMG shortwave band orders differ (RRTMG band 14 becomes RRTMGP band 1, others shift by one) (`GEOS_SolarGridComp.F90`, module header comment on `band_output_supported`).

## 7. Known traps
1. **Output exists only when requested.** Most exports are filled only if the pointer is associated; `no-aerosol` calculations run only if one of the no-aerosol exports is requested (`Run`). A new export that needs new work needs the same guard.
2. **Two code paths for the same number.** Heartbeat quantities come from `UPDATE_EXPORT` (normalized flux times insolation); refresh quantities come from `SORADCORE`. A new flux export needs an internal (normalized) field and a line in `UPDATE_EXPORT`.
3. **Load balancing packs by import.** `SORADCORE` loops over all import names in the import state, so a new import is packed automatically; `AERO` is the only non-field import allowed (`_ASSERT`), a field without horizontal dimensions (`PREF`) is skipped, a horizontal-only field is one slice and a horizontal-plus-vertical field is assumed 3-D (`SORADCORE`, `INPUT_VARS_1`, `SolarFieldGetDims`).
4. **Internals are packed too**, with special cases by name for the aerosol-dependent fields; a new internal must fit the cases in `INT_VARS_1` (`SORADCORE`).
5. **Night points** take the `FILL` value of the internal row (`SORADCORE`, unpack; `Solar_StateSpecs.rc`; a comment in the file says the code must be kept in step with the `FILL` column).
6. **The per-band aerosol arrays are sized by `NUM_BANDS_SOLAR`**, not by what the provider has tables for (parent document, section 6).
7. **RRTMGP needs the NRLSSI2 file.** `SORADCORE` asserts `SOLAR_CYCLE_FILE_NAME` /= `/dev/null` and `USE_NRLSSI2` for RRTMGP.
8. **`ENABLE_SOLAR_RADVAL` changes the restart contents** (module header comment).
9. **Restarts are seamless only at refresh times** (module header comment).
10. **`CALLED_LAST` is read in `Run` and also controls when the alarm should ring**; a source comment states the alarm is always set for the last-step case, so `CALLED_LAST: 0` is not consistent with it (`Run`, comment at "Determine calling sequence").

## 8. Decisions
### Normalized fluxes plus per-step scaling
- Pattern: a full transfer is done at the refresh with the mean cosine zenith over the next interval; fluxes are stored divided by insolation; each step multiplies by the current insolation (`SORADCORE`, `UPDATE_EXPORT`; header comment). [Verified]
- Likely reason: cost; the full calculation is expensive, and insolation changes fastest. [Inferred]
- Question for maintainers: the header says the update "lacks" the path-length effect of the updated sun position and all changes in albedo and atmosphere between refreshes; is that the intended accuracy at the default `SOLAR_DT`?

### Internals that are also exports
- Pattern: refresh-time results are internal rows with `ADD2EXPORT` (`FRIENDLYTO` on the MAPL2 line) (`SetServices`, comment; `Solar_StateSpecs.rc`). [Verified]
- Likely reason: replay corrector reloading (the code comment says so). [Verified from comment]
- Question for maintainers: the comment notes that exports not requested still take internal storage; is the fixed internal state size, independent of the requested exports, intended?

### One file for three schemes
- Pattern: Chou-Suarez, RRTMG and RRTMGP branches are in `SORADCORE` and `UPDATE_EXPORT` of one source file (`SCHEME`). [Verified]
- Likely reason: shared packing, load balancing and export logic. [Inferred]

### Provider-supplied aerosol optics
- Pattern: see the parent document, section 6 and decisions. [Verified]

## 9. Observations (noted, not fixed) and open questions
1. The module header lists the key `AVGR` and describes `RUN_DT`, `DT`; the code reads `SOLAR_DT` (`Initialize`) and no longer reads `SOLARAvrg` or `RUN_DT`.
2. The header "BUGS" list says aerosol properties come from a global method `Get_AeroOptProp` keyed by field name; the code asks the provider through the `AERO` state (`Run`).
3. `NB_CHOU`, `NB_RRTMG`, `NB_RRTMGP` are declared at module level and again in `Run`; a comment says the `Run` copies are redundant.
4. `Run`: the comment on `CALLED_LAST` calls the setting a "kludge" and says the alarm is only correct for the last-step case.
5. The `NUM_BANDS` mismatch message prints a header "SORAD RRTMG:" for all three flags (`Run`).
6. The relative humidity given to the aerosol provider is `QV / MAPL_EQSAT(T, PL)` (`Run`); `IRRAD` uses `GEOS_QSAT`. Not compared.
7. `UPDATE_EXPORT` can load the RRTMGP gas tables itself if the first call precedes any refresh (`UPDATE_EXPORT`).
8. Open question for maintainers: what is the role of `SUN_FLAG` (passed as `TIME` to `MAPL_SunGetInsolation`; default 0)? Which values are in use?
9. In `UPDATE_EXPORT` the comment above the declaration of `rrtmgp_state_set` says the flag is "never set to .true." (issue `GEOSradiation_GridComp#110`) and that the band limits are reloaded on every call; the code sets `rrtmgp_state_set = .true.` in two places in the same routine. *Unverified:* whether the stale-cache problem the comment describes remains.
