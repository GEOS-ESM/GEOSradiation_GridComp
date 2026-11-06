# GEOSradiation_GridComp Architecture

**Audience:** coding agents and developers working on the `GEOSradiation_GridComp` repository (the `RADIATION` child of Physics: solar, infrared and satellite-simulator components).
**Describes:** the MAPL3 line, branch `release/MAPL-v3`.
**Verified against:** `GEOSradiation_GridComp` @ `13f8eb8` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `f50622d`).
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

> **Build status on this line:** built, and run in standalone regression cases.
> On `release/MAPL-v3`, `CMakeLists.txt` builds the `RADIATION` composite
> (`GEOS_RadiationGridComp.F90`) as a shared library with the subdirectories
> `GEOS_RadiationShared`, `GEOSsolar_GridComp`, `GEOSirrad_GridComp` and
> `GEOSsatsim_GridComp`, and adds `regression/`. All four children are ported
> to MAPL3 (`mapl3-porting-notes.md` at the root and in each child directory
> record the porting steps; this document does not repeat them). The parent
> `GEOSphysics_GridComp/CMakeLists.txt` still carries the comment "restricted
> to IRRAD only" on its `GEOSradiation_GridComp` entry; that comment predates
> the ports (*Inferred*), and the nested `CMakeLists.txt` is what decides what is built. The
> Physics library itself is not compiled on this line (the `esma_add_library` block in `../CMakeLists.txt` is commented out, with the comment "Temporarily bypass GEOSphysics top-level library"), so
> `GEOS_PhysicsGridComp.F90`, which still holds the MAPL2 code that adds
> `RADIATION` and connects it to its siblings, is not part of the build; section 7
> describes it as it is in the file.

Paths are relative to this repository's root unless they start with `src/`
(fixture-relative). Statements without a label are Verified and carry a
`(path, routine-or-key)` citation. *Inferred:* and *Unverified:* mark the weaker cases.
Parent documents: Physics (`../../ARCHITECTURE.md` in the fixture) (run order, couplers, synchronization),
GEOSgcm_GridComp (`../../../../docs/ARCHITECTURE.md` in the fixture), fixture (`../../../../../../../docs/ARCHITECTURE.md` in the fixture).
Children: [`GEOSsolar_GridComp/ARCHITECTURE.md`](../GEOSsolar_GridComp/ARCHITECTURE.md),
[`GEOSirrad_GridComp/ARCHITECTURE.md`](../GEOSirrad_GridComp/ARCHITECTURE.md).
Aerosol provider side: Chemistry (`../../@GEOSchem_GridComp/docs/ARCHITECTURE.md` in the fixture), GOCART (`../../@GEOSchem_GridComp/@GOCART/docs/ARCHITECTURE.md` in the fixture).
Where `AGCM.rc` comes from: GEOSgcm_App (`../../../../../../Applications/@GEOSgcm_App/docs/ARCHITECTURE.md` in the fixture).

## 1. Placement
- Fixture key `GEOSradiation_GridComp`, `local: ./src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/@GEOSradiation_GridComp`, `branch: release/MAPL-v3`, `develop: develop` (`components.yaml`, `GEOSradiation_GridComp`).
- Parent: the MAPL2 code in the Physics file adds it as child `RADIATION` with `MAPL_AddChild(GC, NAME='RADIATION', SS=RadiationSetServices, ...)` (`GEOSagcm_GridComp/GEOSphysics_GridComp/GEOS_PhysicsGridComp.F90`, `SetServices`); that file is not compiled on this line (build status above). The ported library exports a routine `Radiation_SetServices` (bottom of `GEOS_RadiationGridComp.F90`) that the standalone regression cases use as the `radiation_setservices_` entry of a configurable Root (`regression/radiation/root.yaml`). The directory is excluded by the `GEOSgcm_GridComp` `.gitignore` (`GEOSphysics_GridComp/.gitignore`, `/@GEOSradiation_GridComp`).
- Nested repository: `RRTMGP` (rte-rrtmgp) at `GEOS_RadiationShared/@RRTMGP`, tag `geos/v1.8+1.0.0`, sparse file `config/RRTMGP.sparse` (`components.yaml`, `RRTMGP`). Its sources are listed by this repository's CMake wrapper, not built from its own build files; the wrapper adds `${MKL_LIBRARIES}` to its dependencies on this line (`GEOS_RadiationShared/RRTMGP_cmake/CMakeLists.txt`).
- Code owners: the radiation team; `CMakeLists.txt` files also the CMake team (`.github/CODEOWNERS`).

## 2. Directory map
| Path | What it is | Module | Document |
|---|---|---|---|
| `GEOS_RadiationGridComp.F90`, `Radiation_StateSpecs.rc` | `RADIATION` composite: adds children, combines their fluxes into heating rates; ACG input for its own state | `GEOS_RadiationGridCompMod` | section 3 |
| `GEOSsolar_GridComp/` | Shortwave (solar) child `SOLAR`; schemes Chou-Suarez `sorad.F90`, RRTMG (`RRTMG/rrtmg_sw`), RRTMGP; ACG file `Solar_StateSpecs.rc` | `GEOS_SolarGridCompMod` | [`GEOSsolar_GridComp/ARCHITECTURE.md`](../GEOSsolar_GridComp/ARCHITECTURE.md) |
| `GEOSirrad_GridComp/` | Longwave (infrared) child `IRRAD`; schemes Chou-Suarez `irrad.F90`, RRTMG (`RRTMG/rrtmg_lw`), RRTMGP; ACG file `Irrad_StateSpecs.rc` | `GEOS_IrradGridCompMod` | [`GEOSirrad_GridComp/ARCHITECTURE.md`](../GEOSirrad_GridComp/ARCHITECTURE.md) |
| `GEOSsatsim_GridComp/` | Optional satellite-simulator child `SATSIM` (COSP: ISCCP, MODIS, MISR, radar, lidar); ACG file `Satsim_StateSpecs.rc` | `GEOS_SatsimGridCompMod` | section 8 |
| `GEOS_RadiationShared/` | Code used by solar and infrared: scheme choice, cloud sampling, constants; builds the `RRTMGP` library | `rad_utils`, `cloud_subcol_gen`, `cloud_condensate_inhomogeneity`, `gettau`, `rad_constants`, `rad_types` | section 9 |
| `bands.txt` | Fortran-comment text listing the Chou-Suarez solar and infrared bands; no source or CMake reference to it was found | | |
| `regression/` | Standalone regression cases `irrad-sa`, `solar-sa`, `satsim-sa`, `radiation` and the `FakeGOCART` aerosol stand-in | | section 10, `regression/README.md` |
| `mapl3-porting-notes.md` (root and in `GEOSirrad_GridComp/`, `GEOSsolar_GridComp/`, `GEOSsatsim_GridComp/`) | Porting notes of the MAPL3 ports: API findings, step lists, baseline setup. Read them for history and rationale; they are not kept in step with the code | | |
| `.circleci/`, `.github/` | CI, code owners, label enforcement | | section 10 |

The child modules are named in the `use ..., only: ...SetServices` lines of the composite (`GEOS_RadiationGridComp.F90`, module header, `use GEOS_IrradGridCompMod, only: irrad_setservices => SetServices` and the two like it). The child directories are listed in `alldirs` (`CMakeLists.txt`).

## 3. The `RADIATION` composite
Source: `GEOS_RadiationGridComp.F90` (about 360 lines). It has no internal state. Its own import and export specs are in the ACG file `Radiation_StateSpecs.rc`, generated by `mapl_acg` in `CMakeLists.txt` (`IMPORT_SPECS EXPORT_SPECS GET_POINTERS DECLARE_POINTERS`); `SetServices` includes only the generated `Radiation_Import___.h` and `Radiation_Export___.h`, and `Run` fetches its pointers by hand with `MAPL_StateGetPointer` (`GEOS_RadiationGridComp.F90`, `SetServices`, `Run`).

### 3.1 Children and how they are added
`SetServices` adds children with `MAPL_GridCompAddChild(gc, name, setservices, config_file, ...)`; each child gets a configuration file of the given name (`GEOS_RadiationGridComp.F90`, `SetServices`):

| Name | Config file | When |
|---|---|---|
| `IRRAD` | `irrad.yaml` | always |
| `SOLAR` | `solar.yaml` | always |
| `SATSIM` | `satsim.yaml` | if `USE_SATSIM + USE_SATSIM_ISCCP + USE_SATSIM_MODIS + USE_SATSIM_RADAR + USE_SATSIM_LIDAR + USE_SATSIM_MISR > 0` (each read with default 0; `get_use_satsim`) |

The regression cases supply these files (`regression/radiation/irrad.yaml`, `solar.yaml`, `satsim.yaml`); where a full-model run gets them was not determined here (Unverified; the Physics code that would add `RADIATION` is not compiled on this line). The MAPL2 `gcm_setup` and `AGCM.rc.tmpl` are not used by the regression cases, which set `USE_SATSIM_ISCCP: 1` and the other `USE_SATSIM*` keys to 0 in `regression/radiation/radiation.yaml`.

### 3.2 Entry points
- `SetServices`: registers `Initialize` and `Run` (`MAPL_GridCompSetEntryPoint`, `Run` with `phase_name="run"`); adds the children; makes the connections and re-exports below; includes the generated import and export specs (`GEOS_RadiationGridComp.F90`, `SetServices`). There is one run phase and no Finalize (MAPL default). The routine has no `MAPL_GenericSetServices` call.
  - `MAPL_GridCompAddConnection` from `IRRAD` to `"<self>"` of `FLX, FLC, FLXA, FLA, DSFDTS0, SFCEM0, TSREFF`, and from `SOLAR` to `"<self>"` of `FSW, FSC, FSWNA, FSCNA`: these fill the composite's own imports (declared in `Radiation_StateSpecs.rc`) from the children's exports, so `Run` can combine them. A connection also makes the child produce the field, which replaces the MAPL2 `MAPL_GetPointer(GEX(child), ..., alloc=.TRUE.)` pattern (`mapl3-porting-notes.md`, "Parent Run reading a child's export state").
  - `MAPL_GridCompAddConnection` of `MCOSZ` from `SOLAR` to `SATSIM`, only if `SATSIM` exists.
  - `MAPL_GridCompReexport` of the exports listed in 3.3.
- `Initialize`: a custom entry point; it does not call a generic initialize (MAPL3 does that itself). It sets two module-level settings that the infrared and solar children both read: `set_inhomogeneity(ih)` from `RAD_CONDENSATE_INHOMOGENEITY` (default 1: 0 homogeneous, 1 beta, 2 gamma) and `initialize_cloud_subcol_gen(...)` with the correlation-length keys `ADL_AM1`, `ADL_AM2`, `ADL_AM30`, `ADL_AM4`, `RDL_AM1`, `RDL_AM2`, `RDL_AM30`, `RDL_AM4` (defaults from `cloud_subcol_gen`), all read with `MAPL_GridCompGetResource` (`GEOS_RadiationGridComp.F90`, `Initialize`; `GEOS_RadiationShared/cloud_subcol_gen.F90`, `def_aam1` etc.). The `IRRAD` and `SOLAR` `Initialize` routines also call `set_inhomogeneity` themselves, with a comment that this is for standalone runs and should eventually be removed (`GEOS_IrradGridComp.F90`, `Initialize`; `GEOS_SolarGridComp.F90`, `Initialize`).
- `Run` (one phase), in order (`GEOS_RadiationGridComp.F90`, `Run`):
  1. Stop the `TOTAL` timer, run the children one at a time with `MAPL_GridCompRunChild`: `IRRAD`, then `SOLAR`, then `SATSIM` if `get_use_satsim` is positive; restart the timer. SATSIM runs after SOLAR because it imports `MCOSZ`. (The MAPL2 code ran all children with `MAPL_GenericRunChildren` in child-index order, `SOLAR`, `IRRAD`, `SATSIM`.)
  2. Get pointers to the imports (`PLEINST` and the connected child exports) and exports. The Edge (interface) imports have Fortran bounds `1:LM+1` from MAPL; `Run` re-bounds them to `0:LM` through a scratch pointer (`p3d`) before using them (`mapl3-porting-notes.md`, "Edge (VLOC=E) import pointers").
  3. Combine: `BLW = DSFDTS0`, `ALW = SFCEM0 - DSFDTS0*TSREFF`, `RADSRF = FSW(LM) + FLX(LM)` (net downward shortwave plus net downward longwave at the surface), `DTDT = ((FLX(0:LM-1) - FLX(1:LM)) + (FSW(0:LM-1) - FSW(1:LM))) * g/cp`. The eight per-species heating rates `RADLW`, `RADSW`, `RADLWC`, `RADSWC`, `RADLWNA`, `RADSWNA`, `RADLWCNA`, `RADSWCNA` are the flux divergence times `g/(cp*dp)` with `dp` from the import `PLEINST` (K s-1). Fluxes are net downward, indexed `0:LM` from the top. Each export is computed only if it is requested (`associated`).

### 3.3 State
- Imports (`Radiation_StateSpecs.rc`): `PLEINST` (instantaneous interface pressure, `Pa`), plus the eleven fields connected from the children (`FLX`, `FLC`, `FLXA`, `FLA`, `DSFDTS0`, `SFCEM0`, `TSREFF`, `FSW`, `FSC`, `FSWNA`, `FSCNA`; the file says they are not meant to be filled by an external coupler). All have `RESTART = SKIP`. All other inputs are the children's imports, which a parent receives as its own imports (the regression cases connect each name once at the `RADIATION` level and it reaches every child that imports it: `regression/README.md`, "`radiation`"; `regression/radiation/root.yaml`).
- Exports declared here (`Radiation_StateSpecs.rc`): `DTDT` (`Pa K s-1`, pressure-weighted), the eight heating rates (`K s-1`), `RADSRF`, `ALW`, `BLW`. Re-exported from `SOLAR` with `MAPL_GridCompReexport`: `DRPAR`, `DFPAR`, `DRNIR`, `DFNIR`, `DRUVR`, `DFUVR` and their `...N` normalized versions, `FCLD`, `TAUCLI`, `TAUCLW`, `CLDTT`, `ALBEDO`, `FSWBAND`, `FSWBANDNA`, and `DROBIO`, `DFOBIO` when `USE_OCEANOBIOGEOCHEM` /= 0 (default 0; read in `SetServices`). Re-exported from `IRRAD`: `LWS`, `LWS0`. All other child exports are not exports of `RADIATION`; `regression/radiation/history.yaml` shows them named under the child component, for example `SATSIM.<name>`.
- Restarts: the composite has none of its own. The children's internal states are checkpointed (the regression cases restart and checkpoint the `IRRAD` and `SOLAR` internals: `regression/README.md`; `regression/radiation/irrad.yaml`); see the child documents.

## 4. Scheme selection (RRTMGP, RRTMG, Chou-Suarez)
- Solar and infrared are chosen independently, by four resource keys read as real numbers where non-zero means on: `USE_RRTMGP_SORAD`, `USE_RRTMG_SORAD`, `USE_RRTMGP_IRRAD`, `USE_RRTMG_IRRAD`, read with `MAPL_GridCompGetResource(gc, ...)` (`GEOS_RadiationShared/rad_utils.F90`, `choose_solar_scheme`, `choose_irrad_scheme`; both take the gridded component `gc` on this line).
- Precedence: RRTMGP, else RRTMG, else Chou-Suarez. Chou-Suarez has no key of its own: it is chosen when both keys of that band are absent or zero (same routines; their comment: "Chou-Suarez is the default").
- `SOLAR` calls `choose_solar_scheme` in `SetServices` (the number of bands sets the ungridded dimension of exported per-band fluxes). `IRRAD` does not call its chooser in `SetServices`: there it sets `USE_RRTMGP = .true.`, `USE_RRTMG = .false.`, `USE_CHOU = .false.` with the comment "TODO: determine from configuration; default to RRTMGP for now", before including the generated specs, so the specs that depend on those flags (the band exports and internals `OLRBbbRG`, `DOLRBbbRGDT`, condition `USE_RRTMG .or. USE_RRTMGP` in `Irrad_StateSpecs.rc`) are declared whichever scheme the keys select. Both children call both choosers in `Run`, to check `NUM_BANDS` (`GEOSsolar_GridComp/GEOS_SolarGridComp.F90`, `SetServices`, `Run`; `GEOSirrad_GridComp/GEOS_IrradGridComp.F90`, `SetServices`, `Run`).
- Band counts: solar 8 (Chou-Suarez), 14 (RRTMG), 14 (RRTMGP); infrared 10, 16, 16. `NUM_BANDS` (read with `MAPL_GridCompGetResource`, no default) must equal solar plus infrared bands of the chosen pair, or both children fail with `Total number of radiation bands is inconsistent!` (`GEOS_SolarGridComp.F90`, `Run`; `GEOS_IrradGridComp.F90`, `Run`). The regression cases set `NUM_BANDS: 30` (14 + 16, right for an RRTMG or RRTMGP pair) with both RRTMGP keys on (`regression/radiation/irrad.yaml`, `solar.yaml`); the MAPL2-style template `AGCM.rc.tmpl` in the GEOSgcm_App directory also sets `NUM_BANDS: 30`. A Chou-Suarez pair (8 + 10 = 18) would need a different value.
- The following two bullets describe the GEOSgcm_App scripts as they are in the fixture; they write resource files in the MAPL2 style (keys with a trailing colon), and whether a full MAPL3 run reads them was not checked. `gcm_setup` writes the RRTMG pair when `AGCM_LM == 72` and the RRTMGP pair otherwise, by deleting the other block of `AGCM.rc.tmpl` (`gcm_setup`, `RRTMG_RADIATION`, `RRTMGP_RADIATION`; `AGCM.rc.tmpl`, "Set RADIATION Parameterizations"). For RRTMGP it also adds the data files `RRTMGP_GAS_LW`, `RRTMGP_GAS_SW`, `RRTMGP_CLOUD_OPTICS_LW`, `RRTMGP_CLOUD_OPTICS_SW` and `SOLAR_LB_MAX_PASSES: 10` (`AGCM.rc.tmpl`). The regression cases give the RRTMGP v1.8 data files in `irrad.yaml` and `solar.yaml`.
- `gcm_run.j` rewrites the word `RRTMG` to `RRTMGP` in `*_instance*.rc` when the RRTMGP block is active (`gcm_run.j`, "Move GOCART to use RRTMGP Bands"). *Inferred:* these are the aerosol-optics tables of the GOCART instances, so that their bands match RRTMGP.
- Differences visible to other components: Chou-Suarez infrared does not use a 3-D `CO2` (falls back to -1); RRTMG infrared sets the no-aerosol fluxes to `MAPL2_UNDEF` (`GEOS_IrradGridComp.F90`, `Update_Flx`, `LW_Driver`); band-resolved diagnostics `OSRBbbRG`, `ISRBbbRG`, `TBRBbbRG`, `OLRBbbRG` exist only for RRTMG and RRTMGP. Details in the child documents.

## 5. Radiation refresh versus every-step update
Both children run on every call of the composite. Each does a full radiative transfer only when its own alarm rings (the refresh) and otherwise updates its exports cheaply from the internal state of the last refresh. On this line each child creates its own alarm in its `Initialize` (MAPL3 no longer creates it for a component; the code comments say the set-up copies MAPL2's `handle_clock_and_main_alarm`): `solar_run_alarm` with interval `SOLAR_DT` (the component name plus `_DT`), `irrad_lw_alarm` with `IRRAD_DT`, `satsim_alarm` with `SATSIM_DT`. Each defaults to the clock time step; `Initialize` of `IRRAD` and `SATSIM` asserts the interval is positive and a multiple of the time step. `SOLAR` `Initialize` has no such assertion: it reads `SOLAR_DT` and sets the alarm interval directly (`GEOS_SolarGridComp.F90`, `Initialize`; searched for `_ASSERT` on the interval). The reference time of the alarm comes from `REFERENCE_DATE`, `REFERENCE_TIME` (default midnight of the current day) and `RUN_AT_INTERVAL_START` (default false) (`GEOSsolar_GridComp/GEOS_SolarGridComp.F90`, `Initialize`; `GEOSirrad_GridComp/GEOS_IrradGridComp.F90`, `Initialize`; `GEOSsatsim_GridComp/GEOS_SatsimGridComp.F90`, `Initialize`). `Run` fetches the alarm from the clock by name. `AGCM.rc.tmpl` sets `SOLAR_DT` from `@SOLAR_DT` (a `set_resolution_params.csh` value in the MAPL2 scripts); the regression cases set `IRRAD_DT`, `SOLAR_DT` and `SATSIM_DT` to 1200 s, equal to their single step. The `RADIATION` composite itself has no alarm of its own.

| | Refresh (alarm rings) | Every call |
|---|---|---|
| `SOLAR` | gas, cloud and aerosol optics; fluxes normalized by top-of-atmosphere insolation, for the mean sun position over the next interval, stored in the internal state | sun position for the current clock step; fluxes = stored normalized fluxes times current insolation (`UPDATE_EXPORT`) |
| `IRRAD` | full longwave transfer; net fluxes and their derivative with respect to surface temperature stored in the internal state | fluxes = stored fluxes plus derivative times (instantaneous surface temperature minus the one used at the refresh) (`Update_Flx`) |
| `SATSIM` | everything, when the alarm rings (nothing otherwise) | none |

Evidence: `GEOSsolar_GridComp/GEOS_SolarGridComp.F90`, `Run`, `UPDATE_EXPORT`, `SORADCORE`; `GEOSirrad_GridComp/GEOS_IrradGridComp.F90`, `Run`, `Update_Flx`, `LW_Driver`; `GEOSsatsim_GridComp/GEOS_SatsimGridComp.F90`, `Run`.
- Order inside one call is set by `CALLED_LAST` (default 1): the update comes first, then the refresh (solar: `UPDATE_FIRST = CalledLast /= 0`). With 0 the refresh comes first.
- Consequence: the heating rates `DTDT`, `RADLW`, `RADSW` that Physics receives each step are computed from flux profiles that change between refreshes only through insolation (shortwave) and surface temperature (longwave). They do not respond to new temperature, humidity or clouds until the next refresh (the solar source comment says the update "lacks" the effects of changed albedo and atmospheric properties: `GEOS_SolarGridComp.F90`, `Run`, comment at `MAPL_SunGetInsolation` in `SORADCORE`).
- The MAPL2 line gave the children's imports an averaging interval (`SOLARAvrg`, `IRRADAvrg`) and a refresh interval. On this line the import specs (`Solar_StateSpecs.rc`, `Irrad_StateSpecs.rc`, IMPORT) have no such columns, and no source in this repository reads `SOLARAvrg` or `IRRADAvrg`: a refresh sees the import values of the step on which the alarm rings, and nothing in this repository time-averages them. (`AGCM.rc.tmpl` still sets both keys to 0.)
- `PLEINST` and `TSINST` are the instantaneous versions of `PLE` and `TS` used by the every-step work (`GEOS_RadiationGridComp.F90`, `Run`; `GEOS_IrradGridComp.F90`, `Update_Flx`). In the MAPL2 code of the fixture the Agcm connects `PLE` to `PLEINST` (`GEOSagcm_GridComp/GEOS_AgcmGridComp.F90`, `SetServices`) and Physics connects `TS` of Surface to `TSINST` (`GEOS_PhysicsGridComp.F90`, `SetServices`); neither file is compiled on this line. The regression cases supply both from a data component (`regression/radiation/root.yaml`).
- Restart: refresh results are internal state, so a restart from a time that is not a refresh time continues from the stored refresh (solar header comment: restarts are seamless only at the time of a full calculation). The regression cases restart the `IRRAD` and `SOLAR` internals from checkpoints (`regression/README.md`).

## 6. Aerosol optics from Chemistry
- Connection: the MAPL2 code of Physics connects the state item `AERO` from `CHEM` to `RAD` (`GEOS_PhysicsGridComp.F90`, `SetServices`, `SHORT_NAME = (/'AERO'/)`, `SRC_ID = CHEM`, `DST_ID = RAD`); that file is not compiled on this line. `SOLAR` and `IRRAD` each declare `AERO` as an import by hand, because the ACG file cannot express a state item: `MAPL_GridCompAddSpec(..., short_name='AERO', itemtype=MAPL_STATEITEM_STATE)` (`GEOS_SolarGridComp.F90`, `SetServices`; `GEOS_IrradGridComp.F90`, `SetServices`; comments there). The MAPL3 calls give no restart attribute for it. In the regression cases `AERO` comes from `FakeGOCART_GridComp`, which supplies an empty state with `implements_aerosol_optics_method = .false.`, so aerosol optics are skipped (`regression/FakeGOCART_GridComp/FakeGOCART_GridComp.F90`; `regression/README.md`, "Common setup").
- Provider (in the full model; the Chemistry composite is not compiled on this line, see its document): Chemistry exports the `AERO` state of the component named by `AERO_PROVIDER` (default `GOCART2G`; `none` gives an empty state whose attribute `implements_aerosol_optics_method` is set `.false.` by the parent) (`@GEOSchem_GridComp/GEOS_ChemGridComp.F90`, `SetServices`, `Init`). A provider that supports radiation sets `implements_aerosol_optics_method = .true.` on its `AERO` state and attaches an ESMF method to it. GOCART2G, CARMA and GMI attach the method under the label `run_aerosol_optics`, the label Radiation executes (`@GEOSchem_GridComp/@GOCART/ESMF/GOCART2G_GridComp/GOCART2G_GridCompMod.F90`, `Initialize`; `@CARMA/CARMAchem_GridCompMod.F90`; `@GMI/GMIchem_GridCompMod.F90`). MAM, DNA and the GEOS-Chem coupler set the attribute true and attach a method labelled `aerosol_optics` (`@MAM/MAM_GridCompMod.F90`; `DNA_GridComp/DNA_GridCompMod.F90`; `GEOSCHEMchem_GridComp/@geos-chem/Interfaces/GEOS/geos_aerocoupler.F90`); see Observation 8.
- Protocol used by Radiation (refresh only; the every-step update does not call the provider):
  1. Read `implements_aerosol_optics_method`; if false, no aerosols are used (`NA = 0`, zero optical depths).
  2. Read the attribute names `relative_humidity_for_aerosol_optics` and `air_pressure_for_aerosol_optics`; if non-empty, write relative humidity and interface pressure into the provider fields of those names.
  3. For each band: set attribute `band_for_aerosol_optics`, execute `run_aerosol_optics`, read the field names given by `extinction_in_air_due_to_ambient_aerosol`, `single_scattering_albedo_of_ambient_aerosol`, `asymmetry_parameter_of_ambient_aerosol`. Values are clipped: extinction at least 0, single-scattering albedo at most 0.999999, asymmetry at most 0.999. Per the source comments the second field is extinction times albedo and the third is extinction times albedo times asymmetry.
  (`GEOSsolar_GridComp/GEOS_SolarGridComp.F90` and `GEOSirrad_GridComp/GEOS_IrradGridComp.F90`, module procedure `compute_provider_aerosol_optics` in each, called from `Run` and `LW_Driver`.)
- Band numbering is one sequence: solar bands are `1..NB_solar`; infrared bands are `NB_solar + 1 ...`, where `NB_solar` is the band count of the solar scheme in use (`BANDS_SOLAR_OFFSET = 0` in `GEOS_SolarGridComp.F90`, `Run`; `OFFSET` in `GEOS_IrradGridComp.F90`, `LW_Driver`). This is why `NUM_BANDS` must match 4.
- Without-aerosol diagnostics: `SOLAR` runs a second transfer without aerosols only when a no-aerosol export is requested (`FSWNA`, `FSCNA`, `FSWBANDNA` and the other `...NA` fields); otherwise it zeroes those internals. `IRRAD` computes clean fluxes when its no-aerosol exports are requested or when no provider implements optics (`GEOS_SolarGridComp.F90`, `Run`; `GEOS_IrradGridComp.F90`, `LW_Driver`, `calc_clrnoa`).
- The composite's own heating rates `RADSWNA`, `RADLWNA`, `RADSWCNA`, `RADLWCNA` request those child exports (`GEOS_RadiationGridComp.F90`, `Run`).

## 7. Interfaces to Physics siblings
Connections made in Physics `SetServices` (`GEOSagcm_GridComp/GEOSphysics_GridComp/GEOS_PhysicsGridComp.F90`). That file holds the MAPL2 code and is not compiled on this line (build status above); the tables describe what it contains, and the names they list are the imports and exports of the ported `RADIATION` (section 3.3 for exports) and its children. In the regression cases a data component stands in for the siblings (`regression/radiation/root.yaml`).

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
| `ALW`, `BLW`, `DRPARN`, `DFPARN`, `DRNIRN`, `DFNIRN`, `DRUVRN`, `DFUVRN`; `LWS0` as `LWDNSRF`; `DROBIO`, `DFOBIO` if `USE_OCEANOBIOGEOCHEM` /= 0 | `SURFACE` |
| `TAUCLI`, `TAUCLW`, `CLDTT`, `DFPAR`, `DRPAR`; `TAUCLI`, `TAUCLW` also as `TAUCLI_avg24`, `TAUCLW_avg24` | `CHEMISTRY` |

`RADIATION` is the last child in the Physics run order (after Chemistry phase 2): it uses the temperature of stage 3 and the aerosols of Chemistry phase 2 (Physics document, run order; `GEOS_PhysicsGridComp.F90`, `Run`).

## 8. Satellite simulator (`SATSIM`)
`GEOSsatsim_GridComp/` has no `ARCHITECTURE.md`; this section is its description.
- Purpose: diagnostic ISCCP, MODIS, MISR, radar (Quickbeam) and lidar outputs computed with COSP code from grid-mean cloud fields; it feeds no tendency (the composite's `Run` uses no `SATSIM` export) (`GEOSsatsim_GridComp/GEOS_SatsimGridComp.F90`; `GEOS_RadiationGridComp.F90`, `Run`).
- `SetServices` registers `Initialize` and `Run` (phase `run`) and includes the ACG-generated import and export specs from `Satsim_StateSpecs.rc` (19 imports: `T`, `QV`, `FCLD`, `QL`...`QG`, `RL`...`RG`, `PLE`, `ZLE`, `TS`, `MCOSZ`, `FRLAND`, `FROCEAN`; 222 exports, some with a 7-bin ungridded dimension; `mapl_acg` in `GEOSsatsim_GridComp/CMakeLists.txt`). The `SATORB` bundle import, the optional `SatSim.rc` table `Masked_Exports::` and the masking step of `SIM_DRIVER` of the MAPL2 line are removed on this line; a comment in `SetServices` and `GEOSsatsim_GridComp/mapl3-porting-notes.md` (step 4) say why and where the old code is. No internal state and no private state are declared (searched `GEOS_SatsimGridComp.F90`, `Satsim_StateSpecs.rc`). Exports keep the MAPL2 undefined value, `1.0e15` (`MAPL2_UNDEF`).
- `Initialize` creates the alarm `satsim_alarm` with interval `SATSIM_DT` (default the time step; section 5). `Run` does work only when the alarm rings and calls `SIM_DRIVER`; the `USE_SATSIM_*` keys choose which simulators are on and are read again there. Other keys: `SATSIM_NCOLUMNS`, `SATSIM_POINTS_PER_ITERATION` (`GEOS_SatsimGridComp.F90`, `Initialize`, `Run`, `SIM_DRIVER`).
- Inputs: `MCOSZ` from `SOLAR` through the composite connection; *Inferred:* everything else is received as a `RADIATION` import like the other children's imports (the names coincide with those of `SOLAR`), as the regression cases do (`regression/README.md`, "`satsim-sa`").
- Regression: `satsim-sa` and `radiation` run the ISCCP simulator only; the lidar, MISR and MODIS exports are not requested there, and the ungridded ISCCP collection is left out of History (`regression/README.md`).
- Sub-trees: `quickbeam/` (radar), `actsim/` (lidar), `llnl/`, plus COSP files at the top level and `cmor/` tables. `Met_Office_license.txt` covers the COSP code. *Unverified:* origin of each sub-tree.

## 9. Shared code (`GEOS_RadiationShared/`)
- `rad_utils.F90`: `choose_solar_scheme`, `choose_irrad_scheme`, `Tbr_from_band_flux` (brightness temperature from band flux).
- `cloud_subcol_gen.F90`: stochastic cloud subcolumns (McICA) with correlation lengths, `clearCounts_threeBand` (low/middle/high/total cloud fractions). `cloud_condensate_inhomogeneity.F90`: homogeneous, beta or gamma condensate (state set by `Initialize` of the composite).
- `gettau.F90`: cloud optical depth for the Chou-Suarez code (`getvistau`, `getnirtau`, `getirtau`; bodies in `getvistau.code`, `getnirtau.code`, `getirtau.code`, included by the preprocessor). `radconstants.F90` (module `rad_constants`): cloud-optics coefficients.
- `NRLSSI2/`: Python scripts and a README that prepare the NRLSSI2 solar-variability input file (`NRLSSI2/README.md`).
- `RRTMGP_cmake/`: CMake wrapper for the nested `RRTMGP` plus GEOS-specific random-number code (`rng/`): the Intel MKL generator if MKL is found, otherwise Mersenne Twister (`RRTMGP_cmake/CMakeLists.txt`).

## 10. Build, test, CI
- `CMakeLists.txt`: `esma_add_library` with sources `GEOS_RadiationGridComp.F90`, `SUBCOMPONENTS GEOS_RadiationShared GEOSsolar_GridComp GEOSirrad_GridComp GEOSsatsim_GridComp`, dependencies `MAPL GEOS_Shared ESMF::ESMF`, `TYPE SHARED`; `mapl_acg` for `Radiation_StateSpecs.rc`; `add_subdirectory(regression)`.
- Child dependencies (their `CMakeLists.txt`): solar `MAPL GEOS_Shared GEOS_RadiationShared RRTMG_SW RRTMG_SW_mods RRTMGP`; infrared `GEOS_Shared MAPL GEOS_RadiationShared RRTMG RRTMG_mods RRTMGP`; satsim `GEOS_Shared MAPL GEOS_RadiationShared`; shared `GEOS_Shared MAPL RRTMGP`. On this line solar, infrared and satsim each run `mapl_acg` on their own `*_StateSpecs.rc`.
- Options: `ENABLE_SOLAR_RADVAL` (defines `SOLAR_RADVAL` in solar; extra developer diagnostics that enlarge the restart state), `RRTMGP_SINGLE_PRECISION` (defines `RTE_USE_SP`); `HAVE_MKL` is defined when MKL is found (`GEOSsolar_GridComp/CMakeLists.txt`, `GEOSirrad_GridComp/CMakeLists.txt`, `GEOS_RadiationShared/RRTMGP_cmake/CMakeLists.txt`).
- Regression: `regression/cases.txt` lists `irrad-sa`, `solar-sa`, `satsim-sa`, `radiation`; `regression/CMakeLists.txt` adds `FakeGOCART_GridComp` and registers the cases with `esma_add_regression_tests` (ctest label `REGRESSION`). Each case runs one component or the composite under a configurable Root on C12L181 with 6 tasks for one 20-minute step, with inputs from restarted data components, and compares checkpoints with a stored baseline (`regression/run_case.cmake`, `regression/README.md`). A case is skipped when its data directory is missing, and the comparison is done only for IntelLLVM builds (`run_case.cmake`). Baselines are from a MAPL2 run with `AERO_PROVIDER: none` (`mapl3-porting-notes.md`, "MAPL2 baselines for the standalone tests"). No other unit-test target was found.
- CI: builds the GEOSgcm fixture at `release/MAPL-v3` with the nested repositories at their release branches (`run-mepo-develop: false`) and runs `ctest -L REGRESSION` (`.circleci/config.yml`; `.github/workflows/workflow.yml`); pull requests need a 0-diff label (`.github/workflows/enforce-labels.yml`).

## 11. Configuration that changes structure
Read with `MAPL_GridCompGetResource`, by label. The keys are written here without the trailing colon of the MAPL2 resource-file style; the regression cases give them in YAML (`regression/radiation/*.yaml`), and the MAPL2-style template is `src/Applications/@GEOSgcm_App/AGCM.rc.tmpl`. Each child is given a configuration file by `MAPL_GridCompAddChild`; `regression/radiation/radiation.yaml` sets `USE_SATSIM_*` for the composite and says SATSIM inherits them (Unverified: the lookup order of MAPL3 `MAPL_GridCompGetResource`).

| Key | Effect |
|---|---|
| `USE_RRTMGP_SORAD`, `USE_RRTMG_SORAD`, `USE_RRTMGP_IRRAD`, `USE_RRTMG_IRRAD` | scheme (section 4) |
| `NUM_BANDS` | must equal solar plus infrared bands (section 4) |
| `SOLAR_DT`, `IRRAD_DT`, `SATSIM_DT`; `RUN_AT_INTERVAL_START`, `REFERENCE_DATE`, `REFERENCE_TIME` | refresh intervals and alarm timing (section 5). `SOLARAvrg`, `IRRADAvrg` are not read on this line |
| `CALLED_LAST` | update before (1) or after (0) refresh |
| `USE_SATSIM`, `USE_SATSIM_ISCCP`, `..._MODIS`, `..._RADAR`, `..._LIDAR`, `..._MISR` | adds the `SATSIM` child and selects simulators |
| `AERO_PROVIDER` (Chemistry), `CO2`, `CO2_PROVIDER` | aerosol optics, 3-D `CO2` import (`CO2: -2` requires a provider, else `IRRAD` raises an error: `GEOS_IrradGridComp.F90`, `SetServices`) |
| `USE_OCEANOBIOGEOCHEM` | adds `DROBIO`, `DFOBIO` (33 bands) |
| `RAD_CONDENSATE_INHOMOGENEITY`, `ADL_*`, `RDL_*` | cloud sampling settings (3.2) |
| `RRTMGP_GAS_*`, `RRTMGP_CLOUD_OPTICS_*`, `RRTMGP_CLOUD_OPTICS_TYPE_*` (`LUT` or `PADE`) | RRTMGP data and cloud-optics form |
| `SOLAR_CONSTANT`, `SOLAR_CYCLE_FILE_NAME`, `USE_NRLSSI2`, `ISOLVAR` | solar forcing (solar document) |

## 12. Known traps
1. **Scheme keys are read in several places.** `SOLAR` reads the same keys at `SetServices`, and both children again in `Run`, through the shared routines (`IRRAD` `SetServices` does not read them on this line, section 4); a new scheme must be added to both children, to `rad_utils`, and to the band counts and the `NUM_BANDS` test (section 4).
2. **`NUM_BANDS` is not derived.** Changing a scheme in the configuration without changing `NUM_BANDS` stops the model in the first `Run`; the aerosol provider's optics tables must have the matching bands (section 6).
3. **Fluxes at the heartbeat are not a full calculation** (section 5). A change that should affect heating immediately (new cloud scheme output, new albedo) reaches `DTDT` only at the next refresh, except through insolation and surface temperature.
4. **Child exports exist only if requested.** A new composite quantity from a child export needs an import in `Radiation_StateSpecs.rc` and a `MAPL_GridCompAddConnection` to `"<self>"` in `SetServices` (3.2), or the child leaves it unfilled.
5. **Level index.** Fluxes are on interfaces `0:LM`, `0` at the top of the atmosphere; heating uses top minus bottom. MAPL3 creates Edge fields with bounds `1:LM+1`, so every routine that uses the `0:LM` convention re-bounds its Edge pointers first (`Run` here; the child `Run` and `Update_Flx` routines).
6. **Adding an import to `SOLAR`** must respect its load-balancing code: only 2-D horizontal and 3-D horizontal-plus-vertical imports are accepted (solar document, trap 3).
7. **Restart size depends on a compile option** (`ENABLE_SOLAR_RADVAL`): restart files from builds with and without it are not interchangeable (`GEOSsolar_GridComp/GEOS_SolarGridComp.F90`, header comment).
8. **Children share module-level settings** (condensate inhomogeneity, correlation lengths) set once in the composite `Initialize`; they are not per-child.

## 13. Observations (noted, not fixed)
1. The `Initialize` description of the composite says the module-level cloud-sampling settings "should be done away with when there is a better treatment of these tables" (`GEOS_RadiationGridComp.F90`, `Initialize`).
2. Under RRTMG infrared, `FLXA` and `FLA` are set to `MAPL2_UNDEF` (`1.0e15`, a module constant of the infrared component on this line), so `RADLWNA` and `RADLWCNA` are computed as undefined minus undefined and come out as 0 (`GEOS_RadiationGridComp.F90`, `Run`; `GEOS_IrradGridComp.F90`, `LW_Driver`, `USE_RRTMG` branch). Derived by reading; not run.
3. `TAUCRIT` has default 0.30 in `IRRAD` and 0.10 in `SOLAR` (`GEOS_IrradGridComp.F90`, `Update_Flx`; `GEOS_SolarGridComp.F90`, `UPDATE_EXPORT`).
4. The relative humidity given to the aerosol provider is computed with `MAPL_EQSAT` in `SOLAR` and with `GEOS_QSAT(..., PASCALS=.true.)` in `IRRAD` (`GEOS_SolarGridComp.F90`, `Run`; `GEOS_IrradGridComp.F90`, `LW_Driver`). Whether the two agree was not checked.
5. The `IRRAD` and `SOLAR` `Initialize` routines each call `set_inhomogeneity`, which the composite `Initialize` also calls; the comments say this is for standalone runs and "eventually remove" (`GEOS_IrradGridComp.F90`, `GEOS_SolarGridComp.F90`, `Initialize`).
6. The band-count constants are declared three times in the solar code (module level and in `Run`); a comment says two are redundant (`GEOS_SolarGridComp.F90`, `Run`, "PMN these three are redundant").
7. `RRTMG_SW_INI` and `RRTMG_LW_INI` are called in every refresh; a comment in the infrared code says they need to be called only once per run (`GEOS_SolarGridComp.F90`, `SORADCORE`; `GEOS_IrradGridComp.F90`, `LW_Driver`).
8. MAM, DNA and the GEOS-Chem coupler attach the aerosol-optics method as `aerosol_optics`, while `SOLAR` and `IRRAD` execute `run_aerosol_optics` (section 6). GOCART2G attaches `run_aerosol_optics` on the top-level state and `aerosol_optics` on its species children. What the labels of MAM, DNA and GEOS-Chem are meant to be was not determined.
9. The MAPL2 code of Physics connects `OX` from Chemistry to Radiation (`GEOS_PhysicsGridComp.F90`, `SetServices`; not compiled on this line), but no import named `OX` exists in `RADIATION`, `SOLAR`, `IRRAD` or `SATSIM` (searched `GEOS_RadiationGridComp.F90`, `GEOSsolar_GridComp/`, `GEOSirrad_GridComp/`, `GEOSsatsim_GridComp/`). Radiation imports `O3`.

Questions for maintainers are in [`decisions.md`](decisions.md).
