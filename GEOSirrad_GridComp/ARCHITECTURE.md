# GEOSirrad_GridComp (IRRAD): architecture

**Audience:** coding agents and developers working on the longwave (infrared) radiation gridded component.
**Describes:** the MAPL2 line, branch `develop`.
**Verified against:** `GEOSradiation_GridComp` @ `0b402f6` (2026-11-05).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `0f2e0fb`; repository `GEOSradiation_GridComp`, mounted at `src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSradiation_GridComp`): `.../@GEOSradiation_GridComp/GEOSirrad_GridComp`.
Other fixtures (for example GEOSadas, GEOSldas) may place this repository elsewhere; paths inside this repository do not change.

Paths are relative to `GEOSirrad_GridComp/` unless they start with `GEOS_RadiationShared/` (same repository root) or `src/` (fixture-relative).
Parent: [Radiation](../docs/ARCHITECTURE.md) (scheme choice, refresh versus update, aerosol optics protocol, interfaces to Physics siblings: read it first). Sibling: [Solar](../GEOSsolar_GridComp/ARCHITECTURE.md). Repository index: [../AGENTS.md](../AGENTS.md).

## 1. Purpose
`IRRAD` computes longwave (thermal infrared) fluxes through a cloudy, aerosol-laden atmosphere with one of three schemes: Chou-Suarez `IRRAD` (default), RRTMG, or RRTMGP. The full transfer, including its linearization with respect to surface temperature, runs only when the component's alarm rings; every call refreshes the exported fluxes from the stored result and the latest surface temperature (`GEOS_IrradGridComp.F90`, header comment, `Run`).

## 2. Place in the hierarchy
- Parent: `RADIATION`, which adds it as `IRRAD` with `MAPL_AddChild(GC, NAME='IRRAD', SS=irradSetServices)` (`GEOS_RadiationGridComp.F90`, `SetServices`). It has no children.
- Files: `GEOS_IrradGridComp.F90` (module `GEOS_IrradGridCompMod`, about 3760 lines, the component and the RRTMGP and RRTMG drivers); `irrad.F90` (module `irradmod`, Chou-Suarez `IRRAD`); `irradconstants.F90`; `Irrad_StateSpecs.rc` (ACG input); `RRTMG/rrtmg_lw/gcm_model/` (AER RRTMG longwave: `modules/` data, `src/` code; libraries `RRTMG`, `RRTMG_mods`); `IrradDoc94.pdf`, `IrradDoc03.pdf` (Chou-Suarez reports; the module header cites the 2003 report).
- RRTMGP code is not in this directory: it is the nested `RRTMGP` repository built by `GEOS_RadiationShared/RRTMGP_cmake` (parent document, section 1).

## 3. Entry points
- `SetServices` (`GEOS_IrradGridComp.F90`, `SetServices`):
  - allocates the RRTMGP private state and attaches it with `ESMF_UserCompSetInternalState(GC, 'RRTMGP_state', ...)`;
  - reads `RUN_DT:` (required), `<name>_DT:` (refresh interval, default heartbeat) and `<name>Avrg:` (averaging interval, default refresh interval), with `<name>` = `IRRAD`;
  - calls `choose_irrad_scheme`;
  - includes the ACG-generated `Irrad_Import___.h`, `Irrad_Export___.h`, `Irrad_Internal___.h`;
  - adds by hand the state-item import `AERO` and, when `CO2:` is -2, the import `CO2` (below);
  - when `RATS_DIAGNOSTICS:` is set, adds the diagnostic exports and internals (section 4);
  - `MAPL_GridCompSetEntryPoint` for `Run` only; `MAPL_GenericSetServices`. Initialize and Finalize are the MAPL generic ones.
- `Run`, one phase (`GEOS_IrradGridComp.F90`, `Run`), every call:
  1. Scheme choice for both bands (`choose_irrad_scheme`, `choose_solar_scheme`) and the `NUM_BANDS:` test (parent document, section 4). The test fails the run with `Total number of radiation bands is inconsistent!`.
  2. For RRTMG and RRTMGP, find which of the supported bands have `OLRBbbRG` or `TBRBbbRG` requested (`band_output`).
  3. Read `CALLED_LAST:` (default 1). If non-zero, `Update_Flx` runs first.
  4. If the run alarm rings (`ESMF_AlarmIsRinging`), turn it off and run `LW_Driver` (the refresh).
  5. If `CALLED_LAST == 0`, `Update_Flx` runs after the refresh.
- `LW_Driver` (contained in `Run`), per refresh:
  - gets pointers to the imports; computes layer pressure `PL` and relative humidity `RH = Q / GEOS_QSAT(T, PL)`;
  - copies `FCLD` (optionally made binary by `RADLW_BINARY_CLOUDS:`);
  - sets `CO2_FIXED` from `CO2:`: -2 selects the 3-D import (not for Chou-Suarez, which falls back to -1); -1 gives the time-dependent value from `getco2`;
  - reads the `RATS_DIAGNOSTICS:` list once (`first`);
  - sets the aerosol band `OFFSET` equal to the solar band count of the solar scheme in use;
  - builds the cloud water paths and effective radii (with defaults 36, 14, 50, 50, 50 microns where the import equals `MAPL_UNDEF`);
  - finds the high, middle and low cloud-layer boundaries `LCLDMH`, `LCLDLM` from `PREF` (`PRS_LOW_MID_CLOUDS:` default 70000 Pa, `PRS_MID_HIGH_CLOUDS:` default 40000 Pa);
  - obtains aerosol optics from the `AERO` import with `compute_provider_aerosol_optics` (parent document, section 6);
  - calls the chosen scheme (`SCHEME:` block): `USE_CHOU` calls `IRRAD` from `irradmod`; `USE_RRTMGP` runs `PROCESS_RRTMGP_LW_BLOCK` per block of columns in an OpenMP loop; `USE_RRTMG` flips the vertical to bottom-to-top, calls `RRTMG_LW_INI` and `RRTMG_LW`, and flips back;
  - then net flux = downward + (negated upward) for the four flux families, `SFCEM_INT = -SFCEM_INT`, stores the surface temperature used (`TS_INT = TS`), and fills `TAUIR`, `CLDTMP`, `CLDPRS`, and the reference-time exports `TSREFF`, `SFCEM0`, `DSFDTS0`, `LWS0`.
- `Update_Flx` (contained in `Run`), every call:
  - computes `CLDTT` (total cloud fraction from `FCLD` in the three pressure bands);
  - `DELT = TSINST - TS_INT`; for Chou-Suarez and RRTMGP: each net or upward flux export = stored flux + stored derivative * `DELT`; downward fluxes are not changed. Under RRTMG the no-aerosol exports are `MAPL_UNDEF` (below);
  - `OLR`, `OLC`, ... from the TOA interface (`-(FLX + DFDTS*DELT)`, assuming zero downward flux at the top, as the source comment says); surface exports `SFCEM`, `LWS`, `LCS`, `FLNS`, ...; band `OLRBbbRG` and `TBRBbbRG`;
  - also writes `DSFDTS0`, `SFCEM0`, `TSREFF` (a source comment suggests these are redundant with the refresh-time values and "should be deprecated");
  - processes the `RATS_DIAGNOSTICS:` exports.
- `compute_provider_aerosol_optics` (module procedure): protocol in the parent document, section 6. For each of the `NB_IRRAD` bands it sets `band_for_aerosol_optics = OFFSET + band`, executes `run_aerosol_optics` and accumulates extinction, albedo-weighted and asymmetry-weighted fields (clipped at 0, 0.999999, 0.999). When the provider does not implement optics, `NA = 0` and the aerosol arrays stay zero.
- `compute_lw_aer_optics`: for RRTMGP only, renormalizes the provider fields into `tau`, `ssa`, `g` of a two-stream optical-properties object and clamps them (comment: slightly out-of-range values were seen with Intel flags).

### RRTMGP path
`LW_Driver` loads the gas k-distribution once into the private state (`RRTMGP_GAS_LW:`, default `rrtmgp-gas-lw-g128.nc`), checks it is longwave (`source_is_internal`) and has 16 bands. Cloud optics are loaded from `RRTMGP_CLOUD_OPTICS_LW:` (`rrtmgp-clouds-lw.nc`) with `RRTMGP_CLOUD_OPTICS_TYPE_LW:` (`LUT` default or `PADE`) on every refresh that needs cloud optics. The calculation runs on blocks of `RRTMGP_LW_BLOCKSIZE:` columns (default 4). Which of four flux sets (clear and all-sky, with and without aerosol) are computed depends on which exports are requested and on whether the provider implements optics: `export_clrnoa`, `export_allnoa`, `export_clrsky`, `export_allsky` and the `calc_*` flags (`LW_Driver`). When no provider implements optics, the no-aerosol results are also used for the aerosol exports. *Inferred:* this saves a second transfer calculation per refresh.

## 4. State
- Imports: by ACG from `Irrad_StateSpecs.rc` category IMPORT: `PLE`, `T`, `QV` (alias `Q`), `QL`, `QI`, `QR`, `QS`, `QG`, `RL`, `RI`, `RR`, `RS`, `RG`, `O3`, `CH4`, `N2O`, `CFC11`, `CFC12`, `HCFC22`, `FCLD`, `TS`, `EMIS` (all with `AVERAGING_INTERVAL = ACCUMINT`, `REFRESH_INTERVAL = MY_STEP`, `RESTART = SKIP`); `PREF` and `TSINST` (no intervals); `CO2` (condition `DT == -2.0` in the ACG file; see Observation 7). By hand: state item `AERO` (`DATATYPE = MAPL_StateItem`, `RESTART = MAPL_RestartSkip`). Who provides them: parent document, section 7. `SOLAR` has more imports (albedos) and no `EMIS`, `TSINST`.
- Internal state (ACG, category INTERNAL): refresh results `FLX`, `FLC`, `FLA`, `FLXD`, `FLXU`, `FLCD`, `FLCU`, `FLAD`, `FLAU`, `FLXA`, `FLXAD`, `FLXAU` (interfaces, `W m-2`), derivatives `DFDTS`, `DFDTSC`, `DFDTSNA`, `DFDTSCNA`, `SFCEM`, `TS`, and per-band `OLRBbbRG`, `DOLRBbbRGDT` for bands 05 to 11 and 15 (condition `USE_RRTMG .or. USE_RRTMGP`). In Fortran the internals are used as `FLX_INT`, `TS_INT`, ... (ACG aliases).
- Exports (ACG, category EXPORT): net, downward and upward flux profiles in four families: all-sky (`FLX`, `FLXD`, `FLXU`), all-sky no-aerosol (`FLXA`, `FLXAD`, `FLXAU`), clear-sky (`FLC`, `FLCD`, `FLCU`) and clear-sky no-aerosol (`FLA`, `FLAD`, `FLAU`), TOA (`OLR`, `OLRA`, `OLC`, `OLCC5`, `OLA`), surface (`SFCEM`, `SFCEM0`, `LWS`, `LWS0`, `LWSA`, `LCS`, `LCSC5`, `LAS`, `FLNS`, `FLNSNA`, `FLNSC`, `FLNSA`, `DSFDTS`, `DSFDTS0`, `TSREFF`), cloud diagnostics (`CLDTMP`, `CLDPRS`, `TAUIR`, `CLDTT`), refresh-frequency cloud fractions `CLDTTLW`, `CLDHILW`, `CLDMDLW`, `CLDLOLW` (from RRTMG subcolumn counts; the ACG comment says they must not be used for heartbeat subsampling), and band exports `OLRBbbRG`, `TBRBbbRG`.
- Internals that are also exports: `DFDTS` has `ADD2EXPORT` (the ACG comment says it is thereby also exported under the same name). The `DFDTS_<species>` internals for RATS also use `add2export = .true.`. Unlike `SOLAR`, the other refresh results are separate internals and exports, filled by `Update_Flx` each call. The `DFDTS` export is therefore the internal itself and holds the last refresh value.
- RATS diagnostics: `RATS_DIAGNOSTICS:` is a list of gas names (for example `CO2 O3`). For each name `SetServices` adds exports `dOLR_*`, `dLWS_*`, `dFLNS_*`, `dSFCEM_*`, `NETTRAP_*`, `COLTRAP_*`, `FLX_*`, `DFDTS_*`, and once `CO2_FIXED`, `DELT`; plus internals `FLXU_RAT`, `FLXD_RAT`, `FLX_RAT`, `DFDTS_RAT`, `SFCEM_RAT` with ungridded dimension equal to the list length. The RRTMG branch calls `RRTMG_LW` once per listed gas with that gas set to zero, then once normally; the differences are the diagnostics (`LW_Driver`, `Update_Flx`). The RRTMGP and Chou-Suarez branches contain no corresponding loop.
- Restarts: internal state is written to `irrad_internal_rst` (`@MODELATMIRRAD_INTERNAL_RESTART_FILE` in `src/Applications/@GEOSgcm_App/AGCM.rc.tmpl`); `AERO` is skipped. The RRTMGP lookup tables are not restarted; they load on the first refresh.
- Private state: type `ty_RRTMGP_state` holds `k_dist` and an `initialized` flag (`GEOS_IrradGridComp.F90`, module). `Update_Flx` reads the RRTMGP band wavenumber limits from it only if it is initialized, so `TBRBbbRG` is `MAPL_UNDEF` before the first refresh.

## 5. Configuration that changes structure
| Key | Effect |
|---|---|
| `USE_RRTMGP_IRRAD:`, `USE_RRTMG_IRRAD:` | scheme; none set means Chou-Suarez (`GEOS_RadiationShared/rad_utils.F90`, `choose_irrad_scheme`) |
| `USE_RRTMGP_SORAD:`, `USE_RRTMG_SORAD:`, `NUM_BANDS:` | solar scheme gives the aerosol band `OFFSET`; `NUM_BANDS:` must equal solar plus infrared band count |
| `IRRAD_DT:`, `IRRADAvrg:`, `CALLED_LAST:` | refresh and averaging intervals, call order |
| `CO2:` | -2: 3-D import (needs `CO2_PROVIDER`, else error in `SetServices`); -1: time-dependent value from `getco2`; positive: fixed |
| `RATS_DIAGNOSTICS:` | adds diagnostic exports and internals; extra `RRTMG_LW` calls (RRTMG only) |
| `RRTMGP_GAS_LW:`, `RRTMGP_CLOUD_OPTICS_LW:`, `RRTMGP_CLOUD_OPTICS_TYPE_LW:`, `RRTMGP_ICE_ROUGHNESS_LW:` (2), `RRTMGP_CLOUD_OVERLAP_TYPE_LW:` (`GEN_MAX_RAN_OVERLAP`), `RRTMGP_LW_BLOCKSIZE:` (4), `RRTMGP_LW_N_GAUSS_ANGLES:`, `RRTMGP_LW_USE_2STREAM:` | RRTMGP data, cloud treatment, blocking; the 2-stream solver is asserted off (`lw_solver_2stream() does not currently support Jacobians`) |
| `RRTMG_ICEFLG:` (3), `RRTMG_LIQFLG:` (1), `RRTMGLW_PARTITION_SIZE:` (4), `RRTMGLW_USE_PRECIP_IN_RADIATION:` (true if `LM > 72`, else false) | RRTMG cloud treatment |
| `RADLW_BINARY_CLOUDS:` | force cloud fraction to 0 or 1 |
| `PRS_LOW_MID_CLOUDS:`, `PRS_MID_HIGH_CLOUDS:`, `TAUCRIT:` (0.30) | cloud-layer split; threshold for `CLDTMP`, `CLDPRS` (the value is divided by 2.13 in the code) |
| Cloud-sampling settings set in `RADIATION` `Initialize` | read by this component through module-level state (parent document, section 3.2) |

## 6. Dependencies
- Libraries (`CMakeLists.txt`): `GEOS_Shared`, `MAPL`, `GEOS_RadiationShared`, `RRTMG`, `RRTMG_mods`, `RRTMGP`; subdirectory `RRTMG`; `HAVE_MKL` is defined when MKL is found. `mapl_acg` generates the import, export, internal, pointer-declaration and pointer-fetch headers from `Irrad_StateSpecs.rc`.
- Modules used: `irradmod`, `rrtmg_lw_rad`, `rrtmg_lw_init`, `parrrtm`, `rrlw_wvn`, `rad_utils` (`Tbr_from_band_flux`, `choose_*_scheme`), `gFTL_StringVector`, `GEOS_UtilsMod` (`GEOS_QSAT`), and inside `LW_Driver` the RRTMGP modules and the cloud sampling code of `GEOS_RadiationShared`. `getco2` comes from `GEOS_Shared`.
- Mixed pairs are not rejected by the code (parent document, section 4).
- RRTMG and RRTMGP both have 16 bands; the band limits are similar but not identical (`GEOS_IrradGridComp.F90`, comment on `band_output_supported`).

## 7. Known traps
1. **Output exists only when requested.** Most exports are filled only if the pointer is associated; under RRTMGP the choice of flux sets to compute depends on the requested exports (`LW_Driver`, `export_*`, `calc_*`). A new export that needs new work needs the same guard.
2. **Two code paths for the same number.** Heartbeat quantities come from `Update_Flx` (stored flux plus derivative times surface-temperature change); refresh quantities come from `LW_Driver`. A new flux export needs a stored internal and a line in `Update_Flx`.
3. **Sign conventions.** Upward fluxes are stored negated (GEOS convention: net downward positive); `DFDTS` is the derivative of the negated upward flux; `SFCEM` is positive after the final sign reversal (`LW_Driver`; comments in `Update_Flx`).
4. **Level index.** Interfaces are `0:LM` from the top; RRTMG is called in reverse order with `1:LM+1` from the bottom, and `LW_Driver` flips on both sides.
5. **Scheme differences in the derivative.** Chou-Suarez gives `DFDTSC = 0` (clear-sky fluxes do not follow the surface temperature between refreshes); RRTMG gives no no-aerosol results and sets `DFDTSNA = DFDTS`, `DFDTSCNA = DFDTSC` (`LW_Driver`).
6. **No load balancing.** Unlike `SOLAR`, there is no day/night packing or work redistribution among processors; columns are independent and RRTMGP blocks use an OpenMP dynamic schedule (`LW_Driver`).
7. **Band output is compile-time restricted.** `OLRBbbRG`, `TBRBbbRG` exist for bands 05-11 and 15 only (`band_output_supported`, ACG file). Adding a band needs the module table, internals for `OLRBbbRG` and `DOLRBbbRGDT`, and exports.
8. **`CALLED_LAST` also decides whether the refresh sees data one step older.** The header says the refresh should occur in the first step of the radiation cycle if the component is called before the users of its fluxes, and during the last step if called after. *Unverified:* which lag results for the default `CALLED_LAST: 1` and the Physics order.
9. **`AERO` dependence.** Provider band numbering assumes `OFFSET = solar band count` (parent document, section 6); a mismatch of the provider tables with `NUM_BANDS:` is not detected here.

## 8. Decisions
### Refresh plus linearized update
- Pattern: a full transfer is done at the refresh and stores the net fluxes and their derivative with respect to surface temperature; each call corrects with the instantaneous surface temperature (`LW_Driver`, `Update_Flx`; header comment). [Verified]
- Likely reason: cost; the surface temperature varies fastest. [Inferred]
- Question for maintainers: only the upward flux is corrected; downward fluxes and atmospheric temperatures are those of the refresh. Is this the intended accuracy at the default `IRRAD_DT`?

### One file for three schemes
- Pattern: Chou-Suarez, RRTMG and RRTMGP calls are in the `SCHEME:` block of `LW_Driver`, with the RRTMGP block code in module procedures of the same file. [Verified]
- Likely reason: shared imports, exports and update logic. [Inferred]

### ACG file for the static state, hand-written code for the dynamic parts
- Pattern: imports, exports and internals that do not depend on runtime keys are in `Irrad_StateSpecs.rc`; `AERO` (a state item), the 3-D `CO2` import and the `RATS_DIAGNOSTICS:` family are added in `SetServices` (ACG file comments). [Verified]
- Likely reason: the MAPL2 ACG has no state-item type and cannot make runtime-dependent names (stated in the ACG comments). [Verified from comment]
- Question for maintainers: `CO2` appears both in the ACG file (condition `DT == -2.0`) and in `SetServices` by hand (Observation 7); which of the two is meant to define the import?

### Provider-supplied aerosol optics
- Pattern: see the parent document, section 6 and decisions. [Verified]

## 9. Observations (noted, not fixed)
1. The module header names the switch `CALL_LAST` ("default `"TRUE"`"); the code reads `CALLED_LAST:` with integer default 1 (`SetServices` header; `Run`).
2. `TAUCRIT:` has default 0.30 here and 0.10 in `SOLAR` (parent document, Observation 3).
3. `RH` for the aerosol provider is computed with `GEOS_QSAT` here and `MAPL_EQSAT` in `SOLAR`; not compared (parent document, Observation 4).
4. `RRTMG_LW_INI` is called in every refresh; a comment in the code says it needs to be done once per run (`LW_Driver`, "pmn: only needs to be done once per run").
5. In the RATS block of the RRTMG branch, the zeroing `select case` has the label `'HCFC22_R'` and the restoring one `'HCFC22'`. A listed `HCFC22` would therefore not be zeroed for its `RRTMG_LW` call but would be "restored" afterwards. Derived by reading; not run. The variable `RATNAMES` (8 names) is declared in `Run` and the list is not used.
6. `RATS_DIAGNOSTICS:` adds exports and internals for every scheme, but only the RRTMG branch computes them (section 4). The pointers `FLX_INT_RAT`, `FLXD_INT_RAT`, `FLXU_INT_RAT`, `DFDTS_RAT` and `SFCEM_INT_RAT` of `Run` are fetched inside the RRTMG branch, while the assignment `FLX_INT_RAT = FLXD_INT_RAT + FLXU_INT_RAT` at the end of `LW_Driver` is guarded only by `nRATS > 0` (`LW_Driver`). *Unverified:* what happens with Chou-Suarez or RRTMGP and a non-empty list (not run).
7. The ACG file declares `CO2` conditionally (`DT == -2.0`), but the generated include is placed in `SetServices` before `DT` is reassigned to the value of `CO2:` (it then holds the averaging interval), and `SetServices` adds `CO2` by hand afterwards when `CO2: -2` (`SetServices`). *Unverified:* the generated header was not inspected (it exists only in a build tree); by reading, the ACG condition is not true at that point.
8. The long name of export `CO2_FIXED` is `lol`, and that of `FLXD_RAT` is `upward_longwave_flux_in_air` (`SetServices`).
9. Chou-Suarez 3-D `CO2` is silently disabled (`CO2_FIXED` falls back to -1) when `CO2: -2` (`LW_Driver`).
10. Open question for maintainers: the DSFDTS0, SFCEM0, TSREFF exports are written both at the refresh and in `Update_Flx` (where the comment says they should be deprecated); the `RADIATION` composite uses `DSFDTS0`, `SFCEM0`, `TSREFF` for `ALW`, `BLW`. Which of the two write sites is the intended source on the step of a refresh?
