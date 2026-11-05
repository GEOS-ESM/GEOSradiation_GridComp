# GEOSradiation_GridComp Architecture

**Audience:** coding agents and developers working on the `GEOSradiation_GridComp` repository (the `RADIATION` child of Physics: solar, infrared and satellite-simulator components).
**Describes:** the MAPL2 line, branch `develop`.
**Verified against:** `GEOSradiation_GridComp` @ `0b402f6` (2026-11-05).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `0f2e0fb`).
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

Paths are relative to this repository's root unless they start with `src/`
(fixture-relative). Statements without a label are Verified and carry a
`(path, routine-or-key)` citation. *Inferred:* and *Unverified:* mark the weaker cases.
Parent documents: [Physics](../../ARCHITECTURE.md) (run order, couplers, synchronization),
[GEOSgcm_GridComp](../../../../docs/ARCHITECTURE.md), [fixture](../../../../../../../docs/ARCHITECTURE.md).
Children: [`GEOSsolar_GridComp/ARCHITECTURE.md`](../GEOSsolar_GridComp/ARCHITECTURE.md),
[`GEOSirrad_GridComp/ARCHITECTURE.md`](../GEOSirrad_GridComp/ARCHITECTURE.md).
Aerosol provider side: [Chemistry](../../@GEOSchem_GridComp/docs/ARCHITECTURE.md), [GOCART](../../@GEOSchem_GridComp/@GOCART/docs/ARCHITECTURE.md).
Where `AGCM.rc` comes from: [GEOSgcm_App](../../../../../../Applications/@GEOSgcm_App/docs/ARCHITECTURE.md).

## 1. Placement
- Fixture key `GEOSradiation_GridComp`, `local: ./src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSradiation_GridComp`, tag `v2.1.0`, `develop: develop` (`components.yaml`, `GEOSradiation_GridComp`).
- Parent: Physics adds it as child `RADIATION` with `MAPL_AddChild(GC, NAME='RADIATION', SS=RadiationSetServices, ...)` (`GEOSagcm_GridComp/GEOSphysics_GridComp/GEOS_PhysicsGridComp.F90`, `SetServices`). The directory is excluded by the `GEOSgcm_GridComp` `.gitignore` (`GEOSphysics_GridComp/.gitignore`, `/@GEOSradiation_GridComp`).
- Nested repository: `RRTMGP` (rte-rrtmgp) at `GEOS_RadiationShared/@RRTMGP`, tag `geos/v1.8+1.0.0`, sparse file `config/RRTMGP.sparse` (`components.yaml`, `RRTMGP`). Its sources are listed by this repository's CMake wrapper, not built from its own build files (`GEOS_RadiationShared/RRTMGP_cmake/CMakeLists.txt`).
- Code owners: the radiation team; `CMakeLists.txt` files also the CMake team (`.github/CODEOWNERS`).

## 2. Directory map
| Path | What it is | Module | Document |
|---|---|---|---|
| `GEOS_RadiationGridComp.F90` | `RADIATION` composite: adds children, combines their fluxes into heating rates | `GEOS_RadiationGridCompMod` | section 3 |
| `GEOSsolar_GridComp/` | Shortwave (solar) child `SOLAR`; schemes Chou-Suarez `sorad.F90`, RRTMG (`RRTMG/rrtmg_sw`), RRTMGP | `GEOS_SolarGridCompMod` | [`GEOSsolar_GridComp/ARCHITECTURE.md`](../GEOSsolar_GridComp/ARCHITECTURE.md) |
| `GEOSirrad_GridComp/` | Longwave (infrared) child `IRRAD`; schemes Chou-Suarez `irrad.F90`, RRTMG (`RRTMG/rrtmg_lw`), RRTMGP; ACG file `Irrad_StateSpecs.rc` | `GEOS_IrradGridCompMod` | [`GEOSirrad_GridComp/ARCHITECTURE.md`](../GEOSirrad_GridComp/ARCHITECTURE.md) |
| `GEOSsatsim_GridComp/` | Optional satellite-simulator child `SATSIM` (COSP: ISCCP, MODIS, MISR, radar, lidar) | `GEOS_SatsimGridCompMod` | section 8 |
| `GEOS_RadiationShared/` | Code used by solar and infrared: scheme choice, cloud sampling, constants; builds the `RRTMGP` library | `rad_utils`, `cloud_subcol_gen`, `cloud_condensate_inhomogeneity`, `gettau`, `rad_constants`, `rad_types` | section 9 |
| `bands.txt` | Fortran-comment text listing the Chou-Suarez solar and infrared bands; no source or CMake reference to it was found | | |
| `.circleci/`, `.github/` | CI, code owners, label enforcement | | section 10 |

The child modules are named in the `use ..., only: ...SetServices` lines of the composite (`GEOS_RadiationGridComp.F90`, module header). The child directories are listed in `alldirs` (`CMakeLists.txt`).

## 3. The `RADIATION` composite
Source: `GEOS_RadiationGridComp.F90` (about 850 lines). It has no internal state and no ACG file (`GEOS_RadiationGridComp.F90`, `SetServices` comment; `CMakeLists.txt` has no `mapl_acg`).

### 3.1 Children and how they are added
`SetServices` adds children with `MAPL_AddChild` (`GEOS_RadiationGridComp.F90`, `SetServices`):

| Name | Module variable | When |
|---|---|---|
| `SOLAR` | `SOL` | always |
| `IRRAD` | `IRR` | always |
| `SATSIM` | `STM` | if `USE_SATSIM + USE_SATSIM_ISCCP + USE_SATSIM_MODIS + USE_SATSIM_RADAR + USE_SATSIM_LIDAR + USE_SATSIM_MISR > 0` (each read with default 0) |

`gcm_setup` substitutes `USE_SATSIM_ISCCP: 1` and the other `USE_SATSIM*` keys 0, so `SATSIM` is a child in experiments made with that script (`src/Applications/@GEOSgcm_App/gcm_setup`, `s/@ISCCP_SATSIM/1/g`; `AGCM.rc.tmpl`, `USE_SATSIM*`). `STM` is set only when `SATSIM` is added.

### 3.2 Entry points
- `SetServices`: registers `Initialize` and `Run` (`MAPL_GridCompSetEntryPoint`); adds the children; declares one import and the exports of 3.3; makes one connection (`MCOSZ` from `SOL` to `STM`, only if `SATSIM` exists); calls `MAPL_GenericSetServices` (`GEOS_RadiationGridComp.F90`, `SetServices`). There is one run phase and no Finalize (MAPL generic default).
- `Initialize`: `MAPL_GenericInitialize` first, then sets two module-level settings that the infrared and solar children both read: `set_inhomogeneity(ih)` from `RAD_CONDENSATE_INHOMOGENEITY:` (default 1: 0 homogeneous, 1 beta, 2 gamma) and `initialize_cloud_subcol_gen(...)` with the correlation-length keys `ADL_AM1:`, `ADL_AM2:`, `ADL_AM30:`, `ADL_AM4:`, `RDL_AM1:`, `RDL_AM2:`, `RDL_AM30:`, `RDL_AM4:` (defaults from `cloud_subcol_gen`) (`GEOS_RadiationGridComp.F90`, `Initialize`; `GEOS_RadiationShared/cloud_subcol_gen.F90`, `def_aam1` etc.). A comment in `Initialize` says they are currently used only by RRTMG and RRTMGP.
- `Run` (one phase), in order (`GEOS_RadiationGridComp.F90`, `Run`):
  1. Get pointers to the exports; for each child export the composite needs, call `MAPL_GetPointer(GEX(child), ..., alloc=.TRUE.)` when the matching composite export is requested. Child exports are filled only if allocated, so this is what makes the child compute them.
  2. `MAPL_GenericRunChildren`: runs `SOLAR`, `IRRAD`, `SATSIM` in child-index order, each followed by the couplers whose source it is (`src/Shared/@MAPL/generic/MAPL_Generic.F90`, `MAPL_GenericRunChildren`).
  3. Combine: `BLW = DSFDTS0`, `ALW = SFCEM0 - DSFDTS0*TSREFF`, `RADSRF = FSW(LM) + FLX(LM)` (net downward shortwave plus net downward longwave at the surface), `DTDT = ((FLX(0:LM-1) - FLX(1:LM)) + (FSW(0:LM-1) - FSW(1:LM))) * g/cp`. The eight per-species heating rates `RADLW`, `RADSW`, `RADLWC`, `RADSWC`, `RADLWNA`, `RADSWNA`, `RADLWCNA`, `RADSWCNA` are the flux divergence times `g/(cp*dp)` with `dp` from the import `PLEINST` (K s-1). Fluxes are net downward, indexed `0:LM` from the top.

### 3.3 State
- Import: `PLEINST` (instantaneous interface pressure, `Pa`) only. All other inputs are the children's imports, which MAPL2 promotes to the composite's import state when no sibling satisfies them (`src/Shared/@MAPL/generic/MAPL_Generic.F90`, `MAPL_WireComponent`).
- Exports declared here: `DTDT` (`Pa K s-1`, pressure-weighted), the eight heating rates (`K s-1`), `RADSRF`, `ALW`, `BLW`. Promoted from `SOLAR` with `MAPL_AddExportSpec(..., CHILD_ID=SOL)`: `DRPAR`, `DFPAR`, `DRNIR`, `DFNIR`, `DRUVR`, `DFUVR` and their `...N` normalized versions, `FCLD`, `TAUCLI`, `TAUCLW`, `CLDTT`, `ALBEDO`, `FSWBAND`, `FSWBANDNA`, and `DROBIO`, `DFOBIO` when `USE_OCEANOBIOGEOCHEM:` /= 0. Promoted from `IRRAD`: `LWS`, `LWS0`. All other child exports are not exports of `RADIATION`; History reads them under the child component names (`src/Applications/@GEOSgcm_App/HISTORY.AGCM.rc.tmpl`, collections that name `'SOLAR'`, `'IRRAD'`, `'SATSIM'`).
- Restarts: the composite has none. The children's internal states are checkpointed (restart files `solar_internal_rst`, `irrad_internal_rst`; keys `@MODELATMSOLAR_INTERNAL_RESTART_FILE`, `@MODELATMIRRAD_INTERNAL_RESTART_FILE` in `AGCM.rc.tmpl`); see the child documents.

## 4. Scheme selection (RRTMGP, RRTMG, Chou-Suarez)
- Solar and infrared are chosen independently, by four resource keys read as real numbers where non-zero means on: `USE_RRTMGP_SORAD:`, `USE_RRTMG_SORAD:`, `USE_RRTMGP_IRRAD:`, `USE_RRTMG_IRRAD:` (`GEOS_RadiationShared/rad_utils.F90`, `choose_solar_scheme`, `choose_irrad_scheme`).
- Precedence: RRTMGP, else RRTMG, else Chou-Suarez. Chou-Suarez has no key of its own: it is chosen when both keys of that band are absent or zero (same routines; their comment: "Chou-Suarez is the default").
- Each child calls its own chooser in `SetServices` (the number of bands sets the ungridded dimension of exported per-band fluxes) and both choosers in `Run`, to check `NUM_BANDS:` (`GEOSsolar_GridComp/GEOS_SolarGridComp.F90`, `SetServices`, `Run`; `GEOSirrad_GridComp/GEOS_IrradGridComp.F90`, `Run`).
- Band counts: solar 8 (Chou-Suarez), 14 (RRTMG), 14 (RRTMGP); infrared 10, 16, 16. `NUM_BANDS:` must equal solar plus infrared bands of the chosen pair, or both children fail with `Total number of radiation bands is inconsistent!` (`GEOS_SolarGridComp.F90`, `Run`; `GEOS_IrradGridComp.F90`, `Run`). `AGCM.rc.tmpl` sets `NUM_BANDS: 30` (14 + 16, right for an RRTMG or RRTMGP pair); a Chou-Suarez pair (8 + 10 = 18) would need a different value in `AGCM.rc`.
- `gcm_setup` writes the RRTMG pair when `AGCM_LM == 72` and the RRTMGP pair otherwise, by deleting the other block of `AGCM.rc.tmpl` (`gcm_setup`, `RRTMG_RADIATION`, `RRTMGP_RADIATION`; `AGCM.rc.tmpl`, "Set RADIATION Parameterizations"). For RRTMGP it also adds the data files `RRTMGP_GAS_LW:`, `RRTMGP_GAS_SW:`, `RRTMGP_CLOUD_OPTICS_LW:`, `RRTMGP_CLOUD_OPTICS_SW:` and `SOLAR_LB_MAX_PASSES: 10` (`AGCM.rc.tmpl`).
- `gcm_run.j` rewrites the word `RRTMG` to `RRTMGP` in `*_instance*.rc` when the RRTMGP block is active (`gcm_run.j`, "Move GOCART to use RRTMGP Bands"). *Inferred:* these are the aerosol-optics tables of the GOCART instances, so that their bands match RRTMGP.
- Differences visible to other components: Chou-Suarez infrared does not use a 3-D `CO2` (falls back to the table); RRTMG infrared sets the no-aerosol fluxes to `MAPL_UNDEF` (`GEOS_IrradGridComp.F90`, `Update_Flx`, `LW_Driver`); band-resolved diagnostics `OSRBbbRG`, `ISRBbbRG`, `TBRBbbRG`, `OLRBbbRG` exist only for RRTMG and RRTMGP. Details in the child documents.

## 5. Radiation refresh versus every-step update
Both children run on every call of the composite. Each does a full radiative transfer only when its own alarm rings (the refresh) and otherwise updates its exports cheaply from the internal state of the last refresh. Alarm interval: the component's `DT:` resource, searched with the component name, so `SOLAR_DT:` and `IRRAD_DT:`; `SATSIM_DT:` for the simulator (`GEOSsolar_GridComp/GEOS_SolarGridComp.F90`, `SetServices`, `SOLAR_DT:`; `src/Shared/@MAPL/generic/MAPL_Generic.F90`, `handle_clock_and_main_alarm`; `src/Shared/@MAPL/base/MAPL_Resource.F90`, `get_labels_with_prefix`; `AGCM.rc.tmpl`, `SOLAR_DT:`). The values come from `set_resolution_params.csh` (3600 s or 1800 s depending on resolution). `RUN_DT:` (the heartbeat) must exist. The `RADIATION` composite itself has no alarm of its own.

| | Refresh (alarm rings) | Every call |
|---|---|---|
| `SOLAR` | gas, cloud and aerosol optics; fluxes normalized by top-of-atmosphere insolation, for the mean sun position over the next interval, stored in the internal state | sun position for the current clock step; fluxes = stored normalized fluxes times current insolation (`UPDATE_EXPORT`) |
| `IRRAD` | full longwave transfer; net fluxes and their derivative with respect to surface temperature stored in the internal state | fluxes = stored fluxes plus derivative times (instantaneous surface temperature minus the one used at the refresh) (`Update_Flx`) |
| `SATSIM` | everything, when the alarm rings (nothing otherwise) | none |

Evidence: `GEOSsolar_GridComp/GEOS_SolarGridComp.F90`, `Run`, `UPDATE_EXPORT`, `SORADCORE`; `GEOSirrad_GridComp/GEOS_IrradGridComp.F90`, `Run`, `Update_Flx`, `LW_Driver`; `GEOSsatsim_GridComp/GEOS_SatsimGridComp.F90`, `Run`.
- Order inside one call is set by `CALLED_LAST:` (default 1): the update comes first, then the refresh (solar: `UPDATE_FIRST = CalledLast /= 0`). With 0 the refresh comes first.
- Consequence: the heating rates `DTDT`, `RADLW`, `RADSW` that Physics receives each step are computed from flux profiles that change between refreshes only through insolation (shortwave) and surface temperature (longwave). They do not respond to new temperature, humidity or clouds until the next refresh (the solar source comment says the update "lacks" the effects of changed albedo and atmospheric properties: `GEOS_SolarGridComp.F90`, `Run`, comment at `MAPL_SunGetInsolation` in `SORADCORE`).
- Imports have `AVERAGING_INTERVAL = SOLARAvrg/IRRADAvrg` and `REFRESH_INTERVAL = SOLAR_DT/IRRAD_DT` (`GEOS_SolarGridComp.F90`, `SetServices`; `Irrad_StateSpecs.rc`, IMPORT). `AGCM.rc.tmpl` sets both averaging keys to 0. A MAPL2 import spec comparison skips the interval test when either averaging interval is 0 (`src/Shared/@MAPL/generic/VarSpec.F90`, `MAPL_VarSpecEQ`). *Inferred:* with 0, a Physics-level connection into `RADIATION` is not forced to use an averaging coupler by the intervals alone, so the refresh sees the value at that step; with a non-zero value a coupler accumulates between refreshes.
- `PLEINST` and `TSINST` are the instantaneous versions of `PLE` and `TS` used by the every-step work (`GEOS_RadiationGridComp.F90`, `Run`; `GEOS_IrradGridComp.F90`, `Update_Flx`). The Agcm connects `PLE` to `PLEINST` (`GEOSagcm_GridComp/GEOS_AgcmGridComp.F90`, `SetServices`); Physics connects `TS` of Surface to `TSINST` (`GEOS_PhysicsGridComp.F90`, `SetServices`).
- Restart: refresh results are internal state, so a restart from a time that is not a refresh time continues from the stored refresh (solar header comment: restarts are seamless only at the time of a full calculation).

## 6. Aerosol optics from Chemistry
- Connection: Physics connects the state item `AERO` from `CHEM` to `RAD` (`GEOS_PhysicsGridComp.F90`, `SetServices`, `SHORT_NAME = (/'AERO'/)`, `SRC_ID = CHEM`, `DST_ID = RAD`). `SOLAR` and `IRRAD` each declare `AERO` as an import with `DATATYPE = MAPL_StateItem`, `RESTART = MAPL_RestartSkip` (`GEOS_SolarGridComp.F90`, `SetServices`; `GEOS_IrradGridComp.F90`, `SetServices`).
- Provider: Chemistry exports the `AERO` state of the component named by `AERO_PROVIDER:` (default `GOCART2G`; `none` gives an empty state whose attribute `implements_aerosol_optics_method` is set `.false.` by the parent) (`@GEOSchem_GridComp/GEOS_ChemGridComp.F90`, `SetServices`, `Init`). A provider that supports radiation sets `implements_aerosol_optics_method = .true.` on its `AERO` state and attaches an ESMF method to it. GOCART2G, CARMA and GMI attach the method under the label `run_aerosol_optics`, the label Radiation executes (`@GEOSchem_GridComp/@GOCART/ESMF/GOCART2G_GridComp/GOCART2G_GridCompMod.F90`, `Initialize`; `@CARMA/CARMAchem_GridCompMod.F90`; `@GMI/GMIchem_GridCompMod.F90`). MAM, DNA and the GEOS-Chem coupler set the attribute true and attach a method labelled `aerosol_optics` (`@MAM/MAM_GridCompMod.F90`; `DNA_GridComp/DNA_GridCompMod.F90`; `GEOSCHEMchem_GridComp/@geos-chem/Interfaces/GEOS/geos_aerocoupler.F90`); see Observation 8.
- Protocol used by Radiation (refresh only; the every-step update does not call the provider):
  1. Read `implements_aerosol_optics_method`; if false, no aerosols are used (`NA = 0`, zero optical depths).
  2. Read the attribute names `relative_humidity_for_aerosol_optics` and `air_pressure_for_aerosol_optics`; if non-empty, write relative humidity and interface pressure into the provider fields of those names.
  3. For each band: set attribute `band_for_aerosol_optics`, execute `run_aerosol_optics`, read the field names given by `extinction_in_air_due_to_ambient_aerosol`, `single_scattering_albedo_of_ambient_aerosol`, `asymmetry_parameter_of_ambient_aerosol`. Values are clipped: extinction at least 0, single-scattering albedo at most 0.999999, asymmetry at most 0.999. Per the source comments the second field is extinction times albedo and the third is extinction times albedo times asymmetry.
  (`GEOSsolar_GridComp/GEOS_SolarGridComp.F90`, `Run`, "Get optical properties of radiatively active aerosols"; `GEOSirrad_GridComp/GEOS_IrradGridComp.F90`, `compute_provider_aerosol_optics`.)
- Band numbering is one sequence: solar bands are `1..NB_solar`; infrared bands are `NB_solar + 1 ...`, where `NB_solar` is the band count of the solar scheme in use (`BANDS_SOLAR_OFFSET = 0` in `GEOS_SolarGridComp.F90`, `Run`; `OFFSET` in `GEOS_IrradGridComp.F90`, `LW_Driver`). This is why `NUM_BANDS:` must match 4.
- Without-aerosol diagnostics: `SOLAR` runs a second transfer without aerosols only when a no-aerosol export is requested (`FSWNA`, `FSCNA`, `FSWBANDNA` and the other `...NA` fields); otherwise it zeroes those internals. `IRRAD` computes clean fluxes when its no-aerosol exports are requested or when no provider implements optics (`GEOS_SolarGridComp.F90`, `Run`; `GEOS_IrradGridComp.F90`, `LW_Driver`, `calc_clrnoa`).
- The composite's own heating rates `RADSWNA`, `RADLWNA`, `RADSWCNA`, `RADLWCNA` request those child exports (`GEOS_RadiationGridComp.F90`, `Run`).

## 7. Interfaces to Physics siblings
Connections made in Physics `SetServices` (`GEOSagcm_GridComp/GEOSphysics_GridComp/GEOS_PhysicsGridComp.F90`):

| Into `RADIATION` | Fields |
|---|---|
| from `MOIST` | `QV`, `QL`, `QI`, `QR`, `QS`, `QG`, `FCLD`, `RL`, `RI`, `RR`, `RS`, `RG` |
| from `SURFACE` | `ALBVR`, `ALBVF`, `ALBNR`, `ALBNF`, `EMIS`, `TS`, `FRLAND`, `FROCEAN`; `TS` also as `TSINST` |
| from `CHEMISTRY` | `OX`, `O3`, `CH4`, `N2O`, `CFC11`, `CFC12`, `HCFC22`; `CO2` only if `CO2_PROVIDER` is `GOCART` or `RRG` and `CO2: -2`; `AERO` (section 6) |
| by hand | `T`: terminated when `SYNCTQ >= 1` and filled from Turbulence's `SAFUPDATE` in stage 3 (`MAPL_TerminateImport`, `CHILD = RAD`; `Run`) |
| from Agcm | `PLEINST` (from `PLE`) (`GEOS_AgcmGridComp.F90`, `SetServices`). *Inferred:* the remaining child imports (`PLE`, `PREF`, `ZLE`, `SATORB`) reach the children as Physics imports of the same name |

| Out of `RADIATION` | To |
|---|---|
| `DTDT` | Physics tendency block (`TIR`; `DTDT`, `DTDTRAD`) |
| `RADLW`, `RADLWC` | `TURBULENCE` |
| `RADLW`, `RADSW` | `MOIST` |
| `ALW`, `BLW`, `DRPARN`, `DFPARN`, `DRNIRN`, `DFNIRN`, `DRUVRN`, `DFUVRN`; `LWS0` as `LWDNSRF`; `DROBIO`, `DFOBIO` if `USE_OCEANOBIOGEOCHEM:` /= 0 | `SURFACE` |
| `TAUCLI`, `TAUCLW`, `CLDTT`, `DFPAR`, `DRPAR`; `TAUCLI`, `TAUCLW` also as `TAUCLI_avg24`, `TAUCLW_avg24` | `CHEMISTRY` |

`RADIATION` is the last child in the Physics run order (after Chemistry phase 2): it uses the temperature of stage 3 and the aerosols of Chemistry phase 2 (Physics architecture, sections 3.2 and 7).

## 8. Satellite simulator (`SATSIM`)
- Purpose: diagnostic ISCCP, MODIS, MISR, radar (Quickbeam) and lidar outputs computed with COSP code from grid-mean cloud fields; it feeds no tendency (the composite's `Run` uses no `SATSIM` export) (`GEOSsatsim_GridComp/GEOS_SatsimGridComp.F90`; `GEOS_RadiationGridComp.F90`, `Run`).
- `SetServices` registers `Run`, about 20 imports (`T`, `QV`, `FCLD`, `QL`...`QG`, `RL`...`RG`, `PLE`, `ZLE`, `TS`, `MCOSZ`, `FRLAND`, `FROCEAN`, bundle `SATORB`) and about 220 exports; it reads an optional file `SatSim.rc` (table `Masked_Exports::`) to add masked exports (`GEOS_SatsimGridComp.F90`, `SetServices`). It has a private ESMF user state `SatSim_State`; no ACG file and no `MAPL_AddInternalSpec` call were found.
- Simulators run when the run alarm rings (`SATSIM_DT:`, set equal to `IRRAD_DT` by `AGCM.rc.tmpl`); the `USE_SATSIM_*` keys choose which of the simulators are switched on, read again in `SIM_DRIVER` (`GEOS_SatsimGridComp.F90`, `Run`, `SIM_DRIVER`). Other keys: `SATSIM_NCOLUMNS:`, `SATSIM_POINTS_PER_ITERATION:`, `AGCM_GRIDNAME:`.
- Inputs: `MCOSZ` from `SOLAR` through the composite connection; `SATORB` from the Agcm (`ORB` export, `GEOS_AgcmGridComp.F90`, `SetServices`); *Inferred:* everything else is promoted to the Physics import state like the other children's imports (the names coincide with those of `SOLAR`).
- Sub-trees: `quickbeam/` (radar), `actsim/` (lidar), `llnl/`, plus COSP files at the top level and `cmor/` tables. `Met_Office_license.txt` covers the COSP code. *Unverified:* origin of each sub-tree.

## 9. Shared code (`GEOS_RadiationShared/`)
- `rad_utils.F90`: `choose_solar_scheme`, `choose_irrad_scheme`, `Tbr_from_band_flux` (brightness temperature from band flux).
- `cloud_subcol_gen.F90`: stochastic cloud subcolumns (McICA) with correlation lengths, `clearCounts_threeBand` (low/middle/high/total cloud fractions). `cloud_condensate_inhomogeneity.F90`: homogeneous, beta or gamma condensate (state set by `Initialize` of the composite).
- `gettau.F90`: cloud optical depth for the Chou-Suarez code (`getvistau`, `getnirtau`, `getirtau`; bodies in `getvistau.code`, `getnirtau.code`, `getirtau.code`, included by the preprocessor). `radconstants.F90` (module `rad_constants`): cloud-optics coefficients.
- `NRLSSI2/`: Python scripts and a README that prepare the NRLSSI2 solar-variability input file (`NRLSSI2/README.md`).
- `RRTMGP_cmake/`: CMake wrapper for the nested `RRTMGP` plus GEOS-specific random-number code (`rng/`): the Intel MKL generator if MKL is found, otherwise Mersenne Twister (`RRTMGP_cmake/CMakeLists.txt`).

## 10. Build, test, CI
- `CMakeLists.txt`: `esma_add_library` with sources `GEOS_RadiationGridComp.F90`, `SUBCOMPONENTS GEOS_RadiationShared GEOSsolar_GridComp GEOSirrad_GridComp GEOSsatsim_GridComp`, dependencies `MAPL GEOS_Shared ESMF::ESMF`.
- Child dependencies (their `CMakeLists.txt`): solar `MAPL GEOS_Shared GEOS_RadiationShared RRTMG_SW RRTMG_SW_mods RRTMGP`; infrared `GEOS_Shared MAPL GEOS_RadiationShared RRTMG RRTMG_mods RRTMGP`; satsim `GEOS_Shared MAPL GEOS_RadiationShared`; shared `GEOS_Shared MAPL RRTMGP`. Only the infrared child runs `mapl_acg`.
- Options: `ENABLE_SOLAR_RADVAL` (defines `SOLAR_RADVAL` in solar; extra developer diagnostics that enlarge the restart state), `RRTMGP_SINGLE_PRECISION` (defines `RTE_USE_SP`); `HAVE_MKL` is defined when MKL is found (`GEOSsolar_GridComp/CMakeLists.txt`, `GEOSirrad_GridComp/CMakeLists.txt`, `GEOS_RadiationShared/RRTMGP_cmake/CMakeLists.txt`).
- No unit-test target was found in this repository. CI builds the GEOSgcm fixture with `develop` of the nested repositories and runs short AMIP and coupled runs (`.circleci/config.yml`; `.github/workflows/workflow.yml`); pull requests need a 0-diff label (`.github/workflows/enforce-labels.yml`).

## 11. Configuration that changes structure
Read with `MAPL_GetResource` from `AGCM.rc` (template `src/Applications/@GEOSgcm_App/AGCM.rc.tmpl`).

| Key | Effect |
|---|---|
| `USE_RRTMGP_SORAD:`, `USE_RRTMG_SORAD:`, `USE_RRTMGP_IRRAD:`, `USE_RRTMG_IRRAD:` | scheme (section 4) |
| `NUM_BANDS:` | must equal solar plus infrared bands (section 4) |
| `SOLAR_DT:`, `IRRAD_DT:`, `SATSIM_DT:`, `SOLARAvrg:`, `IRRADAvrg:` | refresh and averaging intervals (section 5) |
| `CALLED_LAST:` | update before (1) or after (0) refresh |
| `USE_SATSIM:`, `USE_SATSIM_ISCCP:`, `..._MODIS:`, `..._RADAR:`, `..._LIDAR:`, `..._MISR:` | adds the `SATSIM` child and selects simulators |
| `AERO_PROVIDER:` (Chemistry), `CO2:`, `CO2_PROVIDER:` | aerosol optics, 3-D `CO2` import (`CO2: -2` requires a provider, else `IRRAD` raises an error: `GEOS_IrradGridComp.F90`, `SetServices`) |
| `USE_OCEANOBIOGEOCHEM:` | adds `DROBIO`, `DFOBIO` (33 bands) |
| `RAD_CONDENSATE_INHOMOGENEITY:`, `ADL_*`, `RDL_*` | cloud sampling settings (3.2) |
| `RRTMGP_GAS_*:`, `RRTMGP_CLOUD_OPTICS_*:`, `RRTMGP_CLOUD_OPTICS_TYPE_*:` (`LUT` or `PADE`) | RRTMGP data and cloud-optics form |
| `SOLAR_CONSTANT:`, `SOLAR_CYCLE_FILE_NAME:`, `USE_NRLSSI2:`, `ISOLVAR:` | solar forcing (solar document) |

## 12. Known traps
1. **Scheme keys are read in several places.** Each child reads the same keys at `SetServices` and again in `Run` through the shared routines; a new scheme must be added to both children, to `rad_utils`, and to the band counts and the `NUM_BANDS` test (section 4).
2. **`NUM_BANDS:` is not derived.** Changing a scheme in `AGCM.rc` without changing `NUM_BANDS:` stops the model in the first `Run`; the aerosol provider's optics tables must have the matching bands (section 6).
3. **Fluxes at the heartbeat are not a full calculation** (section 5). A change that should affect heating immediately (new cloud scheme output, new albedo) reaches `DTDT` only at the next refresh, except through insolation and surface temperature.
4. **Child exports exist only if requested.** A new composite quantity from a child export needs `MAPL_GetPointer(..., alloc=.TRUE.)` in the composite `Run` (3.2), or the child leaves it unfilled.
5. **Level index.** Fluxes are on interfaces `0:LM`, `0` at the top of the atmosphere; heating uses top minus bottom.
6. **Adding an import to `SOLAR`** must respect its load-balancing code: only 2-D horizontal and 3-D horizontal-plus-vertical imports are accepted (solar document, section 7).
7. **Restart size depends on a compile option** (`ENABLE_SOLAR_RADVAL`): restart files from builds with and without it are not interchangeable (`GEOSsolar_GridComp/GEOS_SolarGridComp.F90`, header comment).
8. **Children share module-level settings** (condensate inhomogeneity, correlation lengths) set once in the composite `Initialize`; they are not per-child.

## 13. Observations (noted, not fixed)
1. The `Initialize` description says it deals with "Aerosols" and "aerosol optical tables in a global module"; the code sets the cloud-sampling settings (`GEOS_RadiationGridComp.F90`, `Initialize`). The `Run` description says a single aerosol bundle is passed to both children; the code passes the `AERO` state through the Physics connection (`Run` description; section 6).
2. Under RRTMG infrared, `FLXA` and `FLA` are `MAPL_UNDEF`, so `RADLWNA` and `RADLWCNA` are computed as undefined minus undefined and come out as 0 (`GEOS_RadiationGridComp.F90`, `Run`; `GEOS_IrradGridComp.F90`, `Update_Flx`, `USE_RRTMG` branch). Derived by reading; not run.
3. `TAUCRIT:` has default 0.30 in `IRRAD` and 0.10 in `SOLAR` (`GEOS_IrradGridComp.F90`, `Update_Flx`; `GEOS_SolarGridComp.F90`, `UPDATE_EXPORT`).
4. The relative humidity given to the aerosol provider is computed with `MAPL_EQSAT` in `SOLAR` and with `GEOS_QSAT` in `IRRAD` (`GEOS_SolarGridComp.F90`, `Run`; `GEOS_IrradGridComp.F90`, `LW_Driver`). Whether the two agree was not checked.
5. The module header of the simulator describes Moist processes (`GEOS_SatsimGridComp.F90`, header). The infrared header names the key `CALL_LAST`; the code reads `CALLED_LAST:`.
6. The band-count constants are declared three times in the solar code (module level and in `Run`); a comment says two are redundant (`GEOS_SolarGridComp.F90`, `Run`, "PMN these three are redundant").
7. `RRTMG_SW_INI` and `RRTMG_LW_INI` are called in every refresh; a comment in the infrared code says they need to be called only once per run (`GEOS_SolarGridComp.F90`, `SORADCORE`; `GEOS_IrradGridComp.F90`, `LW_Driver`).
8. MAM, DNA and the GEOS-Chem coupler attach the aerosol-optics method as `aerosol_optics`, while `SOLAR` and `IRRAD` execute `run_aerosol_optics` (section 6). GOCART2G attaches `run_aerosol_optics` on the top-level state and `aerosol_optics` on its species children. What the labels of MAM, DNA and GEOS-Chem are meant to be was not determined.
9. Physics connects `OX` from Chemistry to Radiation (`GEOS_PhysicsGridComp.F90`, `SetServices`), but no import named `OX` exists in `RADIATION`, `SOLAR`, `IRRAD` or `SATSIM` (searched `GEOS_RadiationGridComp.F90`, `GEOSsolar_GridComp/`, `GEOSirrad_GridComp/`, `GEOSsatsim_GridComp/`). Radiation imports `O3`.

Questions for maintainers are in [`decisions.md`](decisions.md).
