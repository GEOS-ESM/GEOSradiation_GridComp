#include "MAPL_Generic.h"
#undef RTTOV

! GEOS_SatsimGridCompMod drives the COSP satellite simulators (ISCCP, MODIS,
! MISR, CloudSat radar, CALIPSO lidar) from grid-mean cloud and
! thermodynamic fields imported from the model.
module GEOS_SatsimGridCompMod

#define USE_MAPL_UNDEF

   use ESMF
   use MAPL
   use GEOS_UtilsMod

   use gettau

   use MOD_COSP_TYPES
   use MOD_COSP, only: COSP
   !  use MOD_COSP_Modis_Simulator, only: cosp_modis
   use MOD_COSP_Modis_Simulator

   implicit none
   private

   public SetServices

   ! Private State
   type SatSim_State
      private

      integer :: nmask_vars ! number of masked variables
      character(len=ESMF_MAXSTR), pointer :: export_name(:) => null()
      character(len=ESMF_MAXSTR), pointer :: mask_name(:) => null()
      character(len=ESMF_MAXSTR), pointer :: newvar_name(:) => null()
      logical, pointer :: newvar(:) => null()

   end type SatSim_State

   ! Hook for ESMF
   ! -------------

   type SatSim_Wrap
      type(SatSim_State), pointer :: PTR => null()
   end type SatSim_Wrap

contains

   subroutine SetServices(GC, RC)
      type(ESMF_GridComp), intent(inout) :: GC ! gridded component
      integer, optional :: RC ! return code

      integer :: status

      ! Local derived type aliases

      type(MAPL_MetaComp), pointer :: STATE
      type(ESMF_Config) :: CF

      integer :: DIMS(3)

      !   Local derived type aliases
      !   --------------------------
      type(SatSim_State), pointer :: self => null() ! internal, that is
      type(SatSim_Wrap) :: wrap

      integer :: nLines, nCols, m, vindex
      character(len=ESMF_MAXSTR) :: tmpname
      character(len=ESMF_MAXSTR) :: long_name
      character(len=ESMF_MAXSTR) :: units
      integer :: mapl_dims, vlocation
      integer, pointer :: ungridded_dims(:) => null()
      character(len=ESMF_MAXSTR) :: exportName(1)
      type(MAPL_VarSpec), pointer :: ExportSpec(:) => null()
      logical :: found

      !   Wrap internal state for storing in GC; rename legacyState
      !   -------------------------------------
      allocate(self, _STAT)
      wrap%PTR => self

      ! Set the Run entry point
      ! -----------------------

      call MAPL_GridCompSetEntryPoint(GC, ESMF_METHOD_RUN, Run, &
           _RC)

      ! Get the configuration from the component
      !-----------------------------------------

      call ESMF_GridCompGet(GC, CONFIG=CF, _RC)

      ! !IMPORT STATE:

      call MAPL_AddImportSpec(GC, &
           SHORT_NAME='T', &
           long_name='air_temperature', &
           units='K', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddImportSpec(GC, &
           SHORT_NAME='QV', &
           long_name='specific_humidity', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddImportSpec(GC, &
           SHORT_NAME='FCLD', &
           long_name='cloud_area_fraction', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddImportSpec(GC, &
           SHORT_NAME='RL', &
           long_name='liquid_cloud_particle_effective_radius', &
           units='m', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddImportSpec(GC, &
           SHORT_NAME='RI', &
           long_name='ice_phase_cloud_particle_effective_radius', &
           units='m', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddImportSpec(GC, &
           SHORT_NAME='RG', &
           long_name='graupel_particle_effective_radius', &
           units='m', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddImportSpec(GC, &
           SHORT_NAME='RS', &
           long_name='snow_particle_effective_radius', &
           units='m', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddImportSpec(GC, &
           SHORT_NAME='RR', &
           long_name='rain_cloud_particle_effective_radius', &
           units='m', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddImportSpec(GC, &
           SHORT_NAME='PLE', &
           long_name='Edge_pressures', &
           units='Pa', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationEdge, _RC)

      call MAPL_AddImportSpec(GC, &
           SHORT_NAME='ZLE', &
           long_name='Edge_heights', &
           units='m', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationEdge, _RC)

      call MAPL_AddImportSpec(GC, &
           SHORT_NAME='QL', &
           long_name='mass_fraction_of_cloud_liquid_water', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddImportSpec(GC, &
           SHORT_NAME='QI', &
           long_name='mass_fraction_of_cloud_ice_water', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddImportSpec(GC, &
           SHORT_NAME='QR', &
           long_name='mass_fraction_of_falling_rain', &
           units='kg kg-1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddImportSpec(GC, &
           SHORT_NAME='QS', &
           long_name='mass_fraction_of_falling_snow', &
           units='kg kg-1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddImportSpec(GC, &
           SHORT_NAME='QG', &
           long_name='mass_fraction_of_falling_graupel', &
           units='kg kg-1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddImportSpec(GC, &
           SHORT_NAME='TS', &
           long_name='skin_temperature', &
           units='K', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddImportSpec(GC, &
           SHORT_NAME='MCOSZ', &
           long_name='mean_cosine_of_the_solar_zenith_angle', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddImportSpec(GC, &
           SHORT_NAME='FRLAND', &
           long_name='fraction_of_land', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddImportSpec(GC, &
           SHORT_NAME='FROCEAN', &
           long_name='fraction_of_ocean', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      ! !EXPORT STATE:

      !! Fractions in 49 ISCCP Cloud categories
      !!    9 super categories
      !///////////////////////////////////////////
      !     0.   |    .    |    .    |    .       c
      ! 1        |    .    |    .    |    .       l
      !   180.  .|...VII...|..VIII...|....IX....  o
      ! 2        |    .    |    .    |    .       u
      !   310.  .|.........|.........|..........  d
      ! 3        |    .    |    .    |    .
      !   440. ---------------------------------  t
      ! 4        |    .    |    .    |    .       o
      !   560.  .|....IV...|....V....|....VI....  p
      ! 5        |    .    |    .    |    .
      !   680. ---------------------------------  p
      ! 6        | oa . ob |    .    |    .       r
      !   800.  .|....I....|...II....|...III....  e
      ! 7        | ua . ub |    .    |    .       s
      !   sfc    |    .    |    .    |    .       s
      !      ====================================
      !      0.  x   1.3  3.6  9.4  23.  60.   ->
      !        1   2     3    4    5    6    7
      !            optical depth increasing ->
      !
      !   9 Supercategories:
      !     I=cumulus(Cu), II=stratocumulus(StCu), III=stratus(St)
      !    IV=altocumulus(ACu), V=altostratus(ASt), VI=nimbostratus(NSt)
      !    VII=cirrus(Ci), VIII=cirrostratus(CiSt), IX=Deep convection(Cb)
      !
      !   Supercategories further subvided into high/over(O),{middle(M)}, or
      !   low/under(U) and thin (A) and thick (B)
      !
      !   In addition 7 categories of subvisual clouds (Sub) defined by
      !   P_top with Sub_1 above 180 hPa ... etc
      !
      !   In ISCCP simulator 7x7 frequency for these types is arranged thus:
      !
      !         fq_isccp(  itau(1-7),ipres(1-7) )
      !
      !   So frequencies for 7 subvisible cloud types are:
      !
      !          fq_isccp(  1 , 1 ) = subvis above 180.
      !          fq_isccp(  1 , 2 ) = subvis between 180. and 310.
      !               etc. ..
      !

      !
      ! 4 Cumulus (CU) subcategories
      !    fq_isccp(2:3,6:7)
      !-------------------------------------------------
      !    fq_isccp(2,7)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_CU_UA', &
           long_name='isccp_fraction_of_thin_lower_cumulus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(3,7)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_CU_UB', &
           long_name='isccp_fraction_of_thick_lower_cumulus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(2,6)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_CU_OA', &
           long_name='isccp_fraction_of_thin_higher_cumulus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(3,6)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_CU_OB', &
           long_name='isccp_fraction_of_thick_higher_cumulus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      ! 4 Stratocumulus STCU subcategories
      !    fq_isccp(4:5,6:7)
      !-------------------------------------------------
      !    fq_isccp(4,7)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_STCU_UA', &
           long_name='isccp_fraction_of_thin_lower_stratocumulus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(5,7)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_STCU_UB', &
           long_name='isccp_fraction_of_thick_lower_stratocumulus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(4,6)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_STCU_OA', &
           long_name='isccp_fraction_of_thin_higher_stratocumulus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(5,6)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_STCU_OB', &
           long_name='isccp_fraction_of_thick_higher_stratocumulus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      ! 4 Stratus ST subcategories
      !    fq_isccp(6:7,6:7)
      !-------------------------------------------------
      !    fq_isccp(6,7)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_ST_UA', &
           long_name='isccp_fraction_of_thin_lower_stratus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(7,7)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_ST_UB', &
           long_name='isccp_fraction_of_thick_lower_stratus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(6,6)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_ST_OA', &
           long_name='isccp_fraction_of_thin_higher_stratus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(7,6)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_ST_OB', &
           long_name='isccp_fraction_of_thick_higher_stratus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      ! 4 Altocumulus ACU subcategories
      !    fq_isccp(2:3,4:5)
      !-------------------------------------------------
      !    fq_isccp(2,5)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_ACU_UA', &
           long_name='isccp_fraction_of_thin_lower_altocumulus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(3,5)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_ACU_UB', &
           long_name='isccp_fraction_of_thick_lower_altocumulus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(2,4)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_ACU_OA', &
           long_name='isccp_fraction_of_thin_higher_altocumulus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(3,4)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_ACU_OB', &
           long_name='isccp_fraction_of_thick_higher_altocumulus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      ! 4 Altostratus AST subcategories
      !    fq_isccp(4:5,4:5)
      !-------------------------------------------------
      !    fq_isccp(4,5)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_AST_UA', &
           long_name='isccp_fraction_of_thin_lower_altostratus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(5,5)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_AST_UB', &
           long_name='isccp_fraction_of_thick_lower_altostratus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(4,4)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_AST_OA', &
           long_name='isccp_fraction_of_thin_higher_altostratus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(5,4)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_AST_OB', &
           long_name='isccp_fraction_of_thick_higher_altostratus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      ! 4 Nimbostratus NST subcategories
      !    fq_isccp(6:7,4:5)
      !-------------------------------------------------
      !    fq_isccp(6,5)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_NST_UA', &
           long_name='isccp_fraction_of_thin_lower_nimbostratus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(7,5)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_NST_UB', &
           long_name='isccp_fraction_of_thick_lower_nimbostratus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(6,4)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_NST_OA', &
           long_name='isccp_fraction_of_thin_higher_nimbostratus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(7,4)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_NST_OB', &
           long_name='isccp_fraction_of_thick_higher_nimbostratus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      ! 6 Cirrus Ci subcategories
      !    fq_isccp(2:3,1:3)
      !-------------------------------------------------
      !    fq_isccp(2,3)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_CI_UA', &
           long_name='isccp_fraction_of_thin_lower_cirrus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(3,3)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_CI_UB', &
           long_name='isccp_fraction_of_thick_lower_cirrus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(2,2)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_CI_MA', &
           long_name='isccp_fraction_of_thin_middle_cirrus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(3,2)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_CI_MB', &
           long_name='isccp_fraction_of_thick_middle_cirrus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(2,1)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_CI_OA', &
           long_name='isccp_fraction_of_thin_higher_cirrus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(3,1)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_CI_OB', &
           long_name='isccp_fraction_of_thick_higher_cirrus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      ! 6 Cirrostratus CIST subcategories
      !    fq_isccp(4:5,1:3)
      !-------------------------------------------------
      !    fq_isccp(4,3)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_CIST_UA', &
           long_name='isccp_fraction_of_thin_lower_cirrostratus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(5,3)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_CIST_UB', &
           long_name='isccp_fraction_of_thick_lower_cirrostratus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(4,2)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_CIST_MA', &
           long_name='isccp_fraction_of_thin_middle_cirrostratus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(5,2)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_CIST_MB', &
           long_name='isccp_fraction_of_thick_middle_cirrostratus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(4,1)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_CIST_OA', &
           long_name='isccp_fraction_of_thin_higher_cirrostratus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(5,1)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_CIST_OB', &
           long_name='isccp_fraction_of_thick_higher_cirrostratus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      ! 6 Cumulonimbus CB subcategories
      !    fq_isccp(6:7,1:3)
      !-------------------------------------------------
      !    fq_isccp(6,3)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_CB_UA', &
           long_name='isccp_fraction_of_thin_lower_cumulonimbus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(7,3)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_CB_UB', &
           long_name='isccp_fraction_of_thick_lower_cumulonimbus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(6,2)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_CB_MA', &
           long_name='isccp_fraction_of_thin_middle_cumulonimbus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(7,2)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_CB_MB', &
           long_name='isccp_fraction_of_thick_middle_cumulonimbus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(6,1)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_CB_OA', &
           long_name='isccp_fraction_of_thin_higher_cumulonimbus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(7,1)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_CB_OB', &
           long_name='isccp_fraction_of_thick_higher_cumulonimbus', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      ! 7 subvisible/subdetection (SUBV) subcategories
      !    fq_isccp(1,1:7)
      !-------------------------------------------------
      !    fq_isccp(1,1)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_SUBV1', &
           long_name='isccp_fraction_of_subvisible_cloud_0_180_hPa', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(1,2)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_SUBV2', &
           long_name='isccp_fraction_of_subvisible_cloud_180_310_hPa', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(1,3)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_SUBV3', &
           long_name='isccp_fraction_of_subvisible_cloud_310_440_hPa', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(1,4)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_SUBV4', &
           long_name='isccp_fraction_of_subvisible_cloud_440_560_hPa', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(1,5)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_SUBV5', &
           long_name='isccp_fraction_of_subvisible_cloud_560_680_hPa', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(1,6)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_SUBV6', &
           long_name='isccp_fraction_of_subvisible_cloud_680_800_hPa', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(1,7)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP_SUBV7', &
           long_name='isccp_fraction_of_subvisible_cloud_800_SFC_hPa', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      ! 7 exports fields corresponding to pressure bands, with 7 thickess
      ! classes in each -- this is to partially comply with CFMIP data specs
      !    fq_isccp(:,1:7)
      !-------------------------------------------------
      !    fq_isccp(:,1)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CLISCCP1', &
           long_name='isccp_cloud_area_fraction_0_180_hPa', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           ungridded_dims=(/ 7 /), &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(:,2)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CLISCCP2', &
           long_name='isccp_cloud_area_fraction_180_310_hPa', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           ungridded_dims=(/ 7 /), &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(:,3)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CLISCCP3', &
           long_name='isccp_cloud_area_fraction_310_440_hPa', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           ungridded_dims=(/ 7 /), &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(:,4)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CLISCCP4', &
           long_name='isccp_cloud_area_fraction_440_560_hPa', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           ungridded_dims=(/ 7 /), &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(:,5)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CLISCCP5', &
           long_name='isccp_cloud_area_fraction_560_680_hPa', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           ungridded_dims=(/ 7 /), &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(:,6)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CLISCCP6', &
           long_name='isccp_cloud_area_fraction_680_800_hPa', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           ungridded_dims=(/ 7 /), &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(:,7)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CLISCCP7', &
           long_name='isccp_cloud_area_fraction_800_SFC_hPa', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           ungridded_dims=(/ 7 /), &
           vlocation=MAPL_VLocationNone, _RC)

      ! The next 7 exports (ISCCP[1-7]) are deprecated
      !    fq_isccp(:,1)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP1', &
           long_name='isccp_cloud_area_fraction_0_180_hPa', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           ungridded_dims=(/ 7 /), &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(:,2)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP2', &
           long_name='isccp_cloud_area_fraction_180_310_hPa', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           ungridded_dims=(/ 7 /), &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(:,3)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP3', &
           long_name='isccp_cloud_area_fraction_310_440_hPa', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           ungridded_dims=(/ 7 /), &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(:,4)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP4', &
           long_name='isccp_cloud_area_fraction_440_560_hPa', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           ungridded_dims=(/ 7 /), &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(:,5)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP5', &
           long_name='isccp_cloud_area_fraction_560_680_hPa', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           ungridded_dims=(/ 7 /), &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(:,6)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP6', &
           long_name='isccp_cloud_area_fraction_680_800_hPa', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           ungridded_dims=(/ 7 /), &
           vlocation=MAPL_VLocationNone, _RC)

      !    fq_isccp(:,7)
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ISCCP7', &
           long_name='isccp_cloud_area_fraction_800_SFC_hPa', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           ungridded_dims=(/ 7 /), &
           vlocation=MAPL_VLocationNone, _RC)

      ! Other ISCCP ouputs
      !---------------------------------------------------

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CLTISCCP', &
           long_name='isccp_cloud_area_fraction', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      ! same as CLTISCCP, deprecated in favor of CFMIP short name
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='TCLISCCP', &
           long_name='isccp_cloud_area_fraction', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='PCTISCCP', &
           long_name='isccp_air_pressure_at_cloud_top', &
           units='Pa', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      ! same as PCTISCCP, deprecated in favor of CFMIP short name
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CTPISCCP', &
           long_name='isccp_air_pressure_at_cloud_top', &
           units='Pa', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='ALBISCCP', &
           long_name='isccp_cloud_albedo', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='TBISCCP', &
           long_name='isccp_mean_all_sky_10.5_micron_brightness_temp', &
           units='K', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      ! cloud fraction diagnostics from SCOPS

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='SGFCLD', &
           long_name='summed_subgrid_cloud_fraction_from_scops', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      ! lidar_simulator output - not used by CFMIP
      ! MAPL_VLocationCenter is a guess - afe
      ! ---------------------------------------------

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='LIDARPMOL', &
           long_name='calipso_molecular_attenuated_backscatter_signal_power', &
           units='m-1 sr-1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='LIDARPTOT', &
           long_name='calipso_total_attenuated_backscatter_signal_power', &
           units='m-1 sr-1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='LIDARTAUTOT', &
           long_name='calipso_optical_thickess_integrated_from_top_to_level_z', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='RADARZETOT', &
           long_name='cloudsat_total_reflectivity', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CLCALIPSO', &
           long_name='calipso_total_cloud_fraction', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CLLCALIPSO', &
           long_name='calipso_low_level_cloud_fraction', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CLMCALIPSO', &
           long_name='calipso_mid_level_cloud_fraction', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CLHCALIPSO', &
           long_name='calipso_high_level_cloud_fraction', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CLTCALIPSO', &
           long_name='calipso_total_cloud_fraction', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='PARASOLREFL0', &
           long_name='parasol_reflectance', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           ungridded_dims=(/ PARASOL_NREFL /), &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='PARASOLREFL1', &
           long_name='parasol_reflectance_1', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='PARASOLREFL2', &
           long_name='parasol_reflectance_2', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='PARASOLREFL3', &
           long_name='parasol_reflectance_3', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='PARASOLREFL4', &
           long_name='parasol_reflectance_4', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='PARASOLREFL5', &
           long_name='parasol_reflectance_5', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='RADARLTCC', &
           long_name='cloudsat_calipso_total_cloud_amount', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CLCALIPSO2', &
           long_name='calipso_no_cloudsat_cloud_fraction', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CFADLIDARSR532_01', &
           long_name='calipso_scattering_ratio_cfad_01', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CFADLIDARSR532_02', &
           long_name='calipso_scattering_ratio_cfad_02', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CFADLIDARSR532_03', &
           long_name='calipso_scattering_ratio_cfad_03', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CFADLIDARSR532_04', &
           long_name='calipso_scattering_ratio_cfad_04', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CFADLIDARSR532_05', &
           long_name='calipso_scattering_ratio_cfad_05', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CFADLIDARSR532_06', &
           long_name='calipso_scattering_ratio_cfad_06', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CFADLIDARSR532_07', &
           long_name='calipso_scattering_ratio_cfad_07', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CFADLIDARSR532_08', &
           long_name='calipso_scattering_ratio_cfad_08', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CFADLIDARSR532_09', &
           long_name='calipso_scattering_ratio_cfad_09', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CFADLIDARSR532_10', &
           long_name='calipso_scattering_ratio_cfad_10', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CFADLIDARSR532_11', &
           long_name='calipso_scattering_ratio_cfad_11', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CFADLIDARSR532_12', &
           long_name='calipso_scattering_ratio_cfad_12', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CFADLIDARSR532_13', &
           long_name='calipso_scattering_ratio_cfad_13', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CFADLIDARSR532_14', &
           long_name='calipso_scattering_ratio_cfad_14', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CFADLIDARSR532_15', &
           long_name='calipso_scattering_ratio_cfad_15', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      ! Radar simulator exports
      !

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CLOUDSATCFAD01', &
           long_name='cloudsat_radar_reflectivity_cfad', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CLOUDSATCFAD02', &
           long_name='cloudsat_radar_reflectivity_cfad', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CLOUDSATCFAD03', &
           long_name='cloudsat_radar_reflectivity_cfad', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CLOUDSATCFAD04', &
           long_name='cloudsat_radar_reflectivity_cfad', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CLOUDSATCFAD05', &
           long_name='cloudsat_radar_reflectivity_cfad', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CLOUDSATCFAD06', &
           long_name='cloudsat_radar_reflectivity_cfad', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CLOUDSATCFAD07', &
           long_name='cloudsat_radar_reflectivity_cfad', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CLOUDSATCFAD08', &
           long_name='cloudsat_radar_reflectivity_cfad', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CLOUDSATCFAD09', &
           long_name='cloudsat_radar_reflectivity_cfad', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CLOUDSATCFAD10', &
           long_name='cloudsat_radar_reflectivity_cfad', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CLOUDSATCFAD11', &
           long_name='cloudsat_radar_reflectivity_cfad', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CLOUDSATCFAD12', &
           long_name='cloudsat_radar_reflectivity_cfad', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CLOUDSATCFAD13', &
           long_name='cloudsat_radar_reflectivity_cfad', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CLOUDSATCFAD14', &
           long_name='cloudsat_radar_reflectivity_cfad', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='CLOUDSATCFAD15', &
           long_name='cloudsat_radar_reflectivity_cfad', &
           units='1', &
           DIMS=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, _RC)

      !
      ! MODIS simulator output
      !--------------------------------------------------------------------
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSCLDFRCTTL', &
           long_name='modis_cloud_fraction_total_mean', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSCLDFRCWTR', &
           long_name='modis_cloud_fraction_water_mean', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      ! the following is deprecated, use MDSCLDFRCWTR
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSCLDFRCH2O', &
           long_name='modis_cloud_fraction_water_mean', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSCLDFRCICE', &
           long_name='modis_cloud_fraction_ice_mean', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSCLDFRCHI', &
           long_name='modis_cloud_fraction_high_mean', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSCLDFRCMID', &
           long_name='modis_cloud_fraction_mid_mean', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSCLDFRCLO', &
           long_name='modis_cloud_fraction_low_mean', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSOPTHCKTTL', &
           long_name='modis_optical_thickness_total_mean', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSOPTHCKWTR', &
           long_name='modis_optical_thickness_water_mean', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      ! the following is deprecated, use MDSOPTHCKWTR
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSOPTHCKH2O', &
           long_name='modis_optical_thickness_water_mean', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSOPTHCKICE', &
           long_name='modis_optical_thickness_ice_mean', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSOPTHCKTTLLG', &
           long_name='modis_optical_thickness_total_logmean', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSOPTHCKWTRLG', &
           long_name='modis_optical_thickness_water_logmean', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      ! the following is deprecated, use MDSOPTHCKWTRLG
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSOPTHCKH2OLG', &
           long_name='modis_optical_thickness_water_logmean', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSOPTHCKICELG', &
           long_name='modis_optical_thickness_ice_logmean', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSCLDSZWTR', &
           long_name='modis_cloud_particle_size_water_mean', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      ! the following is deprecated, use MDSCLDSZWTR
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSCLDSZH20', &
           long_name='modis_cloud_particle_size_water_mean', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSCLDSZICE', &
           long_name='modis_cloud_particle_size_ice_mean', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSCLDTOPPS', &
           long_name='modis_cloud_top_pressure_total_mean', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSWTRPATH', &
           long_name='modis_liquid_water_path_mean', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      ! the following is deprecated, use MDSWTRPATH
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSH2OPATH', &
           long_name='modis_liquid_water_path_mean', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSICEPATH', &
           long_name='modis_ice_water_path_mean', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST11', &
           long_name='modis_tau_pressure_histogram_bin_1_1', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST12', &
           long_name='modis_tau_pressure_histogram_bin_1_2', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST13', &
           long_name='modis_tau_pressure_histogram_bin_1_3', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST14', &
           long_name='modis_tau_pressure_histogram_bin_1_4', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST15', &
           long_name='modis_tau_pressure_histogram_bin_1_5', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST16', &
           long_name='modis_tau_pressure_histogram_bin_1_6', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST17', &
           long_name='modis_tau_pressure_histogram_bin_1_7', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST21', &
           long_name='modis_tau_pressure_histogram_bin_2_1', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST22', &
           long_name='modis_tau_pressure_histogram_bin_2_2', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST23', &
           long_name='modis_tau_pressure_histogram_bin_2_3', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST24', &
           long_name='modis_tau_pressure_histogram_bin_2_4', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST25', &
           long_name='modis_tau_pressure_histogram_bin_2_5', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST26', &
           long_name='modis_tau_pressure_histogram_bin_2_6', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST27', &
           long_name='modis_tau_pressure_histogram_bin_2_7', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST31', &
           long_name='modis_tau_pressure_histogram_bin_3_1', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST32', &
           long_name='modis_tau_pressure_histogram_bin_3_2', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST33', &
           long_name='modis_tau_pressure_histogram_bin_3_3', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST34', &
           long_name='modis_tau_pressure_histogram_bin_3_4', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST35', &
           long_name='modis_tau_pressure_histogram_bin_3_5', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST36', &
           long_name='modis_tau_pressure_histogram_bin_3_6', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST37', &
           long_name='modis_tau_pressure_histogram_bin_3_7', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST41', &
           long_name='modis_tau_pressure_histogram_bin_4_1', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST42', &
           long_name='modis_tau_pressure_histogram_bin_4_2', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST43', &
           long_name='modis_tau_pressure_histogram_bin_4_3', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST44', &
           long_name='modis_tau_pressure_histogram_bin_4_4', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST45', &
           long_name='modis_tau_pressure_histogram_bin_4_5', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST46', &
           long_name='modis_tau_pressure_histogram_bin_4_6', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST47', &
           long_name='modis_tau_pressure_histogram_bin_4_7', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST51', &
           long_name='modis_tau_pressure_histogram_bin_5_1', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST52', &
           long_name='modis_tau_pressure_histogram_bin_5_2', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST53', &
           long_name='modis_tau_pressure_histogram_bin_5_3', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST54', &
           long_name='modis_tau_pressure_histogram_bin_5_4', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST55', &
           long_name='modis_tau_pressure_histogram_bin_5_5', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST56', &
           long_name='modis_tau_pressure_histogram_bin_5_6', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST57', &
           long_name='modis_tau_pressure_histogram_bin_5_7', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST61', &
           long_name='modis_tau_pressure_histogram_bin_6_1', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST62', &
           long_name='modis_tau_pressure_histogram_bin_6_2', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST63', &
           long_name='modis_tau_pressure_histogram_bin_6_3', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST64', &
           long_name='modis_tau_pressure_histogram_bin_6_4', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST65', &
           long_name='modis_tau_pressure_histogram_bin_6_5', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST66', &
           long_name='modis_tau_pressure_histogram_bin_6_6', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST67', &
           long_name='modis_tau_pressure_histogram_bin_6_7', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST71', &
           long_name='modis_tau_pressure_histogram_bin_7_1', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST72', &
           long_name='modis_tau_pressure_histogram_bin_7_2', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST73', &
           long_name='modis_tau_pressure_histogram_bin_7_3', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST74', &
           long_name='modis_tau_pressure_histogram_bin_7_4', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST75', &
           long_name='modis_tau_pressure_histogram_bin_7_5', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST76', &
           long_name='modis_tau_pressure_histogram_bin_7_6', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MDSTAUPRSHIST77', &
           long_name='modis_tau_pressure_histogram_bin_7_7', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      !
      ! MISR simulator output
      !--------------------------------------------------------------------
      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MISRMNCLDTP', &
           long_name='MISR_mead_cloud_top_height', &
           units='m', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MISRCLDAREA', &
           long_name='MISR_cloud_area', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MISRLYRTP0', &
           long_name='MISR_layer_top', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MISRLYRTP250', &
           long_name='MISR_layer_top', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MISRLYRTP750', &
           long_name='MISR_layer_top', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MISRLYRTP1250', &
           long_name='MISR_layer_top', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MISRLYRTP1750', &
           long_name='MISR_layer_top', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MISRLYRTP2250', &
           long_name='MISR_layer_top', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MISRLYRTP2750', &
           long_name='MISR_layer_top', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MISRLYRTP3500', &
           long_name='MISR_layer_top', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MISRLYRTP4500', &
           long_name='MISR_layer_top', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MISRLYRTP6000', &
           long_name='MISR_layer_top', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MISRLYRTP8000', &
           long_name='MISR_layer_top', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MISRLYRTP10000', &
           long_name='MISR_layer_top', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MISRLYRTP12000', &
           long_name='MISR_layer_top', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MISRLYRTP14000', &
           long_name='MISR_layer_top', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MISRLYRTP16000', &
           long_name='MISR_layer_top', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MISRLYRTP18000', &
           long_name='MISR_layer_top', &
           units='1', &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MISRFQ0', &
           long_name='MISR_cloud_area', &
           units='1', &
           ungridded_dims=(/ 7 /), &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MISRFQ250', &
           long_name='MISR_cloud_area', &
           units='1', &
           ungridded_dims=(/ 7 /), &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MISRFQ750', &
           long_name='MISR_cloud_area', &
           units='1', &
           ungridded_dims=(/ 7 /), &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MISRFQ1250', &
           long_name='MISR_cloud_area', &
           units='1', &
           ungridded_dims=(/ 7 /), &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MISRFQ1750', &
           long_name='MISR_cloud_area', &
           units='1', &
           ungridded_dims=(/ 7 /), &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MISRFQ2250', &
           long_name='MISR_cloud_area', &
           units='1', &
           ungridded_dims=(/ 7 /), &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MISRFQ2750', &
           long_name='MISR_cloud_area', &
           units='1', &
           ungridded_dims=(/ 7 /), &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MISRFQ3500', &
           long_name='MISR_cloud_area', &
           units='1', &
           ungridded_dims=(/ 7 /), &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MISRFQ4500', &
           long_name='MISR_cloud_area', &
           units='1', &
           ungridded_dims=(/ 7 /), &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MISRFQ6000', &
           long_name='MISR_cloud_area', &
           units='1', &
           ungridded_dims=(/ 7 /), &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MISRFQ8000', &
           long_name='MISR_cloud_area', &
           units='1', &
           ungridded_dims=(/ 7 /), &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MISRFQ10000', &
           long_name='MISR_cloud_area', &
           units='1', &
           ungridded_dims=(/ 7 /), &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MISRFQ12000', &
           long_name='MISR_cloud_area', &
           units='1', &
           ungridded_dims=(/ 7 /), &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MISRFQ14000', &
           long_name='MISR_cloud_area', &
           units='1', &
           ungridded_dims=(/ 7 /), &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MISRFQ16000', &
           long_name='MISR_cloud_area', &
           units='1', &
           ungridded_dims=(/ 7 /), &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddExportSpec(GC, &
           SHORT_NAME='MISRFQ18000', &
           long_name='MISR_cloud_area', &
           units='1', &
           ungridded_dims=(/ 7 /), &
           DIMS=MAPL_DimsHorzOnly, &
           vlocation=MAPL_VLocationNone, _RC)

      call MAPL_AddImportSpec(GC, &
           SHORT_NAME='SATORB', &
           long_name='Satellite_orbits', &
           units='days', &
           DIMS=MAPL_DimsHorzOnly, &
           DATATYPE=MAPL_BundleItem, &
           _RC)

      ! Set the Profiling timers
      ! ------------------------

      call MAPL_TimerAdd(GC, NAME="DRIVER", _RC)
      call MAPL_TimerAdd(GC, NAME="-MISC", _RC)
      call MAPL_TimerAdd(GC, NAME="-COSP", _RC)

      ! parse resource file for masks
      CF = ESMF_ConfigCreate(_RC)
      inquire(file='SatSim.rc', exist=found)
      if (found) then
         call ESMF_ConfigLoadFile(CF, 'SatSim.rc', _RC)
         call ESMF_ConfigGetDim(CF, nLines, nCols, LABEL='Masked_Exports::', RC=status)
         if (status == ESMF_SUCCESS) then
            self%nmask_vars = nLines
            allocate(self%export_name(nLines), self%mask_name(nLines), self%newvar_name(nLines), self%newvar(nLines), &
                 _STAT)
            call ESMF_ConfigFindLabel(CF, 'Masked_Exports::', _RC)
            do m = 1, nLines
               call ESMF_ConfigNextLine(CF, _RC)
               call ESMF_ConfigGetAttribute(CF, self%export_name(m), _RC)
               call ESMF_ConfigGetAttribute(CF, self%mask_name(m), _RC)
               call ESMF_ConfigGetAttribute(CF, tmpname, RC=status)
               if (status == ESMF_SUCCESS) then
                  self%newvar_name(m) = tmpname
                  self%newvar(m) = .true.
               else
                  self%newvar_name(m) = self%export_name(m)
                  self%newvar(m) = .false.
               end if
            end do
         end if
      else
         self%nmask_vars = 0
      end if
      if (MAPL_AM_I_Root()) then
         write(*, *)'Parsing of satsim.rc'
         do m = 1, self%nmask_vars
            write(*, *)m, self%newvar(m)
            write(*, *)trim(self%export_name(m)), ' ', trim(self%mask_name(m)), ' ', trim(self%newvar_name(m))
         end do
      end if

      call MAPL_GetObjectFromGC(GC, STATE, _RC)

      do m = 1, self%nmask_vars
         if (self%newvar(m)) then
            call MAPL_Get(STATE, ExportSpec=ExportSpec, _RC)
            vindex = MAPL_VarSpecGetIndex(ExportSpec, self%export_name(m), _RC)
            call MAPL_VarSpecGet(ExportSpec(vindex), long_name=long_name, &
                 units=units, DIMS=mapl_dims, vlocation=vlocation, ungridded_dims=ungridded_dims, _RC)
            if (associated(ungridded_dims)) then
               call MAPL_AddExportSpec(GC, &
                    SHORT_NAME=self%newvar_name(m), &
                    long_name=long_name, &
                    units=units, &
                    ungridded_dims=ungridded_dims, &
                    DIMS=mapl_dims, &
                    vlocation=vlocation, _RC)
               nullify(ungridded_dims)
            else
               call MAPL_AddExportSpec(GC, &
                    SHORT_NAME=self%newvar_name(m), &
                    long_name=long_name, &
                    units=units, &
                    DIMS=mapl_dims, &
                    vlocation=vlocation, _RC)
            end if
            exportName(1) = self%newvar_name(m)
            call MAPL_DoNotDeferExport(GC, exportName, _RC)
         end if
         call MAPL_DoNotDeferExport(GC, self%export_name, _RC)
      end do

      !   Store internal state in GC
      !   --------------------------
      call ESMF_UserCompSetInternalState(GC, 'SatSim_State', wrap, status)
      _VERIFY(status)

      ! Set generic init and final methods
      ! ----------------------------------

      call MAPL_GenericSetServices(GC, _RC)

      _RETURN(_SUCCESS)

   end subroutine SetServices

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine Run(GC, IMPORT, EXPORT, CLOCK, RC)
      type(ESMF_GridComp), intent(inout) :: GC ! Gridded component
      type(ESMF_State), intent(inout) :: IMPORT ! Import state
      type(ESMF_State), intent(inout) :: EXPORT ! Export state
      type(ESMF_Clock), intent(inout) :: CLOCK ! The clock
      integer, optional, intent(out) :: RC ! Error code:

      integer :: status

      ! Local derived type aliases

      type(MAPL_MetaComp), pointer :: STATE
      type(ESMF_Config) :: CF
      type(ESMF_State) :: INTERNAL
      type(ESMF_Alarm) :: ALARM

      ! Local variables

      integer :: IM, JM, LM
      real, pointer, dimension(:, :) :: LONS
      real, pointer, dimension(:, :) :: LATS

      ! pointers for Export outputs

      ! ISCCP/Icarus

      !CU family
      real, pointer, dimension(:, :) :: ISCCP_CU_OA, ISCCP_CU_OB
      real, pointer, dimension(:, :) :: ISCCP_CU_UA, ISCCP_CU_UB

      !STCU family
      real, pointer, dimension(:, :) :: ISCCP_STCU_OA, ISCCP_STCU_OB
      real, pointer, dimension(:, :) :: ISCCP_STCU_UA, ISCCP_STCU_UB

      !ST family
      real, pointer, dimension(:, :) :: ISCCP_ST_OA, ISCCP_ST_OB
      real, pointer, dimension(:, :) :: ISCCP_ST_UA, ISCCP_ST_UB

      !ACU family
      real, pointer, dimension(:, :) :: ISCCP_ACU_OA, ISCCP_ACU_OB
      real, pointer, dimension(:, :) :: ISCCP_ACU_UA, ISCCP_ACU_UB

      !AST family
      real, pointer, dimension(:, :) :: ISCCP_AST_OA, ISCCP_AST_OB
      real, pointer, dimension(:, :) :: ISCCP_AST_UA, ISCCP_AST_UB

      !NST family
      real, pointer, dimension(:, :) :: ISCCP_NST_OA, ISCCP_NST_OB
      real, pointer, dimension(:, :) :: ISCCP_NST_UA, ISCCP_NST_UB

      !CI family
      real, pointer, dimension(:, :) :: ISCCP_CI_OA, ISCCP_CI_OB
      real, pointer, dimension(:, :) :: ISCCP_CI_MA, ISCCP_CI_MB
      real, pointer, dimension(:, :) :: ISCCP_CI_UA, ISCCP_CI_UB

      !CIST family
      real, pointer, dimension(:, :) :: ISCCP_CIST_OA, ISCCP_CIST_OB
      real, pointer, dimension(:, :) :: ISCCP_CIST_MA, ISCCP_CIST_MB
      real, pointer, dimension(:, :) :: ISCCP_CIST_UA, ISCCP_CIST_UB

      !Cb family
      real, pointer, dimension(:, :) :: ISCCP_CB_OA, ISCCP_CB_OB
      real, pointer, dimension(:, :) :: ISCCP_CB_MA, ISCCP_CB_MB
      real, pointer, dimension(:, :) :: ISCCP_CB_UA, ISCCP_CB_UB

      !SubVisible Family
      real, pointer, dimension(:, :) :: ISCCP_SUBV1, ISCCP_SUBV2
      real, pointer, dimension(:, :) :: ISCCP_SUBV3, ISCCP_SUBV4
      real, pointer, dimension(:, :) :: ISCCP_SUBV5, ISCCP_SUBV6
      real, pointer, dimension(:, :) :: ISCCP_SUBV7

      ! "stacked" fields (im*jm*7)
      real, pointer, dimension(:, :, :) :: CLISCCP1, CLISCCP2
      real, pointer, dimension(:, :, :) :: CLISCCP3, CLISCCP4
      real, pointer, dimension(:, :, :) :: CLISCCP5, CLISCCP6
      real, pointer, dimension(:, :, :) :: CLISCCP7

      ! "stacked" fields (im*jm*7)
      real, pointer, dimension(:, :, :) :: ISCCP1, ISCCP2
      real, pointer, dimension(:, :, :) :: ISCCP3, ISCCP4
      real, pointer, dimension(:, :, :) :: ISCCP5, ISCCP6
      real, pointer, dimension(:, :, :) :: ISCCP7

      !other 2Ds
      real, pointer, dimension(:, :) :: TCLISCCP, CLTISCCP
      real, pointer, dimension(:, :) :: CTPISCCP, PCTISCCP
      real, pointer, dimension(:, :) :: ALBISCCP
      real, pointer, dimension(:, :) :: TBISCCP

      ! LIDAR/Calipso exports

      real, pointer, dimension(:, :) :: RADARLTCC
      real, pointer, dimension(:, :) :: CLLCALIPSO
      real, pointer, dimension(:, :) :: CLMCALIPSO
      real, pointer, dimension(:, :) :: CLHCALIPSO
      real, pointer, dimension(:, :) :: CLTCALIPSO

      real, pointer, dimension(:, :, :) :: PARASOLREFL0
      real, pointer, dimension(:, :) :: PARASOLREFL1
      real, pointer, dimension(:, :) :: PARASOLREFL2
      real, pointer, dimension(:, :) :: PARASOLREFL3
      real, pointer, dimension(:, :) :: PARASOLREFL4
      real, pointer, dimension(:, :) :: PARASOLREFL5

      real, pointer, dimension(:, :, :) :: SGFCLD
      real, pointer, dimension(:, :, :) :: LIDARPMOL, LIDARPTOT
      real, pointer, dimension(:, :, :) :: LIDARTAUTOT
      real, pointer, dimension(:, :, :) :: RADARZETOT
      real, pointer, dimension(:, :, :) :: CLCALIPSO2
      real, pointer, dimension(:, :, :) :: CLCALIPSO

      real, pointer, dimension(:, :, :) :: CFADLIDARSR532_01
      real, pointer, dimension(:, :, :) :: CFADLIDARSR532_02
      real, pointer, dimension(:, :, :) :: CFADLIDARSR532_03
      real, pointer, dimension(:, :, :) :: CFADLIDARSR532_04
      real, pointer, dimension(:, :, :) :: CFADLIDARSR532_05
      real, pointer, dimension(:, :, :) :: CFADLIDARSR532_06
      real, pointer, dimension(:, :, :) :: CFADLIDARSR532_07
      real, pointer, dimension(:, :, :) :: CFADLIDARSR532_08
      real, pointer, dimension(:, :, :) :: CFADLIDARSR532_09
      real, pointer, dimension(:, :, :) :: CFADLIDARSR532_10
      real, pointer, dimension(:, :, :) :: CFADLIDARSR532_11
      real, pointer, dimension(:, :, :) :: CFADLIDARSR532_12
      real, pointer, dimension(:, :, :) :: CFADLIDARSR532_13
      real, pointer, dimension(:, :, :) :: CFADLIDARSR532_14
      real, pointer, dimension(:, :, :) :: CFADLIDARSR532_15

      ! Radar/cloudsat exports

      real, pointer, dimension(:, :, :) :: CLOUDSATCFAD01
      real, pointer, dimension(:, :, :) :: CLOUDSATCFAD02
      real, pointer, dimension(:, :, :) :: CLOUDSATCFAD03
      real, pointer, dimension(:, :, :) :: CLOUDSATCFAD04
      real, pointer, dimension(:, :, :) :: CLOUDSATCFAD05
      real, pointer, dimension(:, :, :) :: CLOUDSATCFAD06
      real, pointer, dimension(:, :, :) :: CLOUDSATCFAD07
      real, pointer, dimension(:, :, :) :: CLOUDSATCFAD08
      real, pointer, dimension(:, :, :) :: CLOUDSATCFAD09
      real, pointer, dimension(:, :, :) :: CLOUDSATCFAD10
      real, pointer, dimension(:, :, :) :: CLOUDSATCFAD11
      real, pointer, dimension(:, :, :) :: CLOUDSATCFAD12
      real, pointer, dimension(:, :, :) :: CLOUDSATCFAD13
      real, pointer, dimension(:, :, :) :: CLOUDSATCFAD14
      real, pointer, dimension(:, :, :) :: CLOUDSATCFAD15

      ! MODIS simulator exports

      real, pointer, dimension(:, :) :: MDSCLDFRCTTL
      real, pointer, dimension(:, :) :: MDSCLDFRCWTR
      real, pointer, dimension(:, :) :: MDSCLDFRCH2O
      real, pointer, dimension(:, :) :: MDSCLDFRCICE
      real, pointer, dimension(:, :) :: MDSCLDFRCHI
      real, pointer, dimension(:, :) :: MDSCLDFRCMID
      real, pointer, dimension(:, :) :: MDSCLDFRCLO
      real, pointer, dimension(:, :) :: MDSOPTHCKTTL
      real, pointer, dimension(:, :) :: MDSOPTHCKWTR
      real, pointer, dimension(:, :) :: MDSOPTHCKH2O
      real, pointer, dimension(:, :) :: MDSOPTHCKICE
      real, pointer, dimension(:, :) :: MDSOPTHCKTTLLG
      real, pointer, dimension(:, :) :: MDSOPTHCKWTRLG
      real, pointer, dimension(:, :) :: MDSOPTHCKH2OLG
      real, pointer, dimension(:, :) :: MDSOPTHCKICELG
      real, pointer, dimension(:, :) :: MDSCLDSZWTR
      real, pointer, dimension(:, :) :: MDSCLDSZH20
      real, pointer, dimension(:, :) :: MDSCLDSZICE
      real, pointer, dimension(:, :) :: MDSCLDTOPPS
      real, pointer, dimension(:, :) :: MDSWTRPATH
      real, pointer, dimension(:, :) :: MDSH2OPATH
      real, pointer, dimension(:, :) :: MDSICEPATH

      real, pointer, dimension(:, :) :: MDSTAUPRSHIST11
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST12
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST13
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST14
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST15
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST16
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST17

      real, pointer, dimension(:, :) :: MDSTAUPRSHIST21
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST22
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST23
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST24
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST25
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST26
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST27

      real, pointer, dimension(:, :) :: MDSTAUPRSHIST31
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST32
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST33
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST34
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST35
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST36
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST37

      real, pointer, dimension(:, :) :: MDSTAUPRSHIST41
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST42
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST43
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST44
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST45
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST46
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST47

      real, pointer, dimension(:, :) :: MDSTAUPRSHIST51
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST52
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST53
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST54
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST55
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST56
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST57

      real, pointer, dimension(:, :) :: MDSTAUPRSHIST61
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST62
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST63
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST64
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST65
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST66
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST67

      real, pointer, dimension(:, :) :: MDSTAUPRSHIST71
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST72
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST73
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST74
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST75
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST76
      real, pointer, dimension(:, :) :: MDSTAUPRSHIST77

      ! MISR simulator exports

      real, pointer, dimension(:, :) :: MISRMNCLDTP
      real, pointer, dimension(:, :) :: MISRCLDAREA
      real, pointer, dimension(:, :) :: MISRLYRTP0
      real, pointer, dimension(:, :) :: MISRLYRTP250
      real, pointer, dimension(:, :) :: MISRLYRTP750
      real, pointer, dimension(:, :) :: MISRLYRTP1250
      real, pointer, dimension(:, :) :: MISRLYRTP1750
      real, pointer, dimension(:, :) :: MISRLYRTP2250
      real, pointer, dimension(:, :) :: MISRLYRTP2750
      real, pointer, dimension(:, :) :: MISRLYRTP3500
      real, pointer, dimension(:, :) :: MISRLYRTP4500
      real, pointer, dimension(:, :) :: MISRLYRTP6000
      real, pointer, dimension(:, :) :: MISRLYRTP8000
      real, pointer, dimension(:, :) :: MISRLYRTP10000
      real, pointer, dimension(:, :) :: MISRLYRTP12000
      real, pointer, dimension(:, :) :: MISRLYRTP14000
      real, pointer, dimension(:, :) :: MISRLYRTP16000
      real, pointer, dimension(:, :) :: MISRLYRTP18000
      real, pointer, dimension(:, :, :) :: MISRFQ0
      real, pointer, dimension(:, :, :) :: MISRFQ250
      real, pointer, dimension(:, :, :) :: MISRFQ750
      real, pointer, dimension(:, :, :) :: MISRFQ1250
      real, pointer, dimension(:, :, :) :: MISRFQ1750
      real, pointer, dimension(:, :, :) :: MISRFQ2250
      real, pointer, dimension(:, :, :) :: MISRFQ2750
      real, pointer, dimension(:, :, :) :: MISRFQ3500
      real, pointer, dimension(:, :, :) :: MISRFQ4500
      real, pointer, dimension(:, :, :) :: MISRFQ6000
      real, pointer, dimension(:, :, :) :: MISRFQ8000
      real, pointer, dimension(:, :, :) :: MISRFQ10000
      real, pointer, dimension(:, :, :) :: MISRFQ12000
      real, pointer, dimension(:, :, :) :: MISRFQ14000
      real, pointer, dimension(:, :, :) :: MISRFQ16000
      real, pointer, dimension(:, :, :) :: MISRFQ18000

      ! Get my internal MAPL_Generic state
      !-----------------------------------

      call MAPL_GetObjectFromGC(GC, STATE, _RC)

      ! Get parameters from generic state.
      !-----------------------------------

      call MAPL_Get(STATE, IM=IM, JM=JM, LM=LM, &
           CF=CF, &
           LONS=LONS, &
           LATS=LATS, &
           INTERNAL_ESMF_STATE=INTERNAL, &
           RUNALARM=ALARM, &
           _RC)

      if (.not. ESMF_AlarmIsRinging(ALARM)) then
         _RETURN(_SUCCESS)
      end if

      call MAPL_TimerOn(STATE, "DRIVER")

      call SIM_DRIVER(IM, JM, LM, _RC)

      call MAPL_TimerOff(STATE, "DRIVER")

      _RETURN(_SUCCESS)

   contains

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine SIM_DRIVER(IM, JM, LM, RC)
         integer, intent(in) :: IM, JM, LM
         integer, optional, intent(out) :: RC

         !  Locals

         character(len=ESMF_MAXSTR) :: NAME
         integer :: status

         ! Import arrays

         ! These are the gen-u-ine import arrays in GEOS-5 dimensions
         ! PLE and ZLE are indexed 0:LM
         real, pointer, dimension(:, :, :) :: T, PLE, QV, FCLD, ZLE
         real, pointer, dimension(:, :, :) :: RDFL, RDFI
         real, pointer, dimension(:, :, :) :: RDFR, RDFS, RDFG
         real, pointer, dimension(:, :, :) :: QLTOT, QITOT, QRTOT, QSTOT, QGTOT
         real, pointer, dimension(:, :) :: MCOSZ, FRLAND, TS, FROCEAN

         ! These are the same, in logical 2 dimensions for icarus's sake
         real, dimension(IM * JM, 0:LM) :: ZLE2D
         real, dimension(IM * JM) :: MCOSZCOSP, FRLANDCOSP, TSCOSP, FROCEANCOSP

         ! These are the same, converted to upside-down and 2d
         real, dimension(IM * JM, LM) :: TCOSP, PLOCOSP, QVCOSP, FCLDCOSP
         real, dimension(IM * JM, LM) :: RDFLCOSP, RDFICOSP
         real, dimension(IM * JM, LM) :: QLTOTCOSP, QITOTCOSP
         real, dimension(IM * JM, LM) :: QRTOTCOSP, QSTOTCOSP, QGTOTCOSP
         real, dimension(IM * JM, LM) :: RDFRCOSP, RDFSCOSP, RDFGCOSP
         real, dimension(IM * JM, 0:LM) :: PLECOSP, PLE2DTOA0INV, ZLECOSP

         ! Local variables needed for all simulators
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         real, dimension(IM, JM, LM) :: PLO, RH, QST
         real, dimension(IM * JM, LM) :: RHCOSP
         integer :: Npoints
         real, dimension(:, :, :), target, allocatable :: frac_out
         !real, dimension(:,:,:),target,allocatable :: frac_outinv
         real, dimension(IM * JM, LM) :: frac_ls
         integer :: isccp_overlap
         logical :: DEBUG_GC
         integer :: i, j, k, l, icb, ict
         integer :: scops_debug = 0
         integer :: idebug = 0
         integer :: idebugcol = 0
         real, dimension(:, :), pointer :: column_frac_out

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !  The follwing assignments are for variables needed for COSP
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         integer :: melt_lay ! melting layer model off=0, on=1
         integer :: Naero ! Number of aerosol species
         integer :: surface_radar ! surface=1,spaceborne=0
         ! Time [days] - need to do something about this - afe
         double precision :: time = 1.D0
         double precision :: time_bands(2, 1) = 1.D0

         ! COSP types

         type(cosp_config) :: cfg ! Configuration options
         type(cosp_gridbox) :: gbx ! Gridbox information. Input for COSP
         type(cosp_subgrid) :: sgx ! Subgrid outputs
         type(cosp_sgradar) :: sgradar ! Output from radar simulator
         type(cosp_sghydro) :: sghydro ! Input to radar simulator
         type(cosp_sglidar) :: sglidar ! Output from lidar simulator
         type(cosp_isccp) :: isccp ! Output from ISCCP simulator
         type(cosp_vgrid) :: vgrid ! Information on vertical grid of stats
         type(cosp_radarstats) :: stradar ! Summary statistics from radar simulator
         type(cosp_lidarstats) :: stlidar ! Summary statistics from lidar simulator
         type(cosp_modis) :: modis ! Output from MODIS simulator
         type(cosp_misr) :: misr ! Output from MISR simulator

         ! RRTOV paramters needed only by COSP, otherwise unused
         ! -------------------------------------------------------

         integer :: Npoints_it
         logical :: use_precipitation_fluxes, use_reff
         integer :: Plat
         integer :: Sat
         integer :: Inst
         integer :: Nchan
         integer :: Ichan(0)
         real :: SurfEm(0)
         real :: ZenAng
         real :: co2, ch4, n2o, co
         !

         !  The follwing definitions are for variables needed for COSP and isccp/icarus
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         real :: CCA(IM, JM, LM) ! convective clouds fraction, goes into scops
         real :: CCACOSP(IM * JM, LM) ! convective clouds fraction, goes into scops
         real :: DTAU_S(IM, JM, LM)
         real :: dtau_sCOSP(IM * JM, LM)
         real, dimension(IM, JM, LM, 4) :: CWC, REFF
         real, dimension(IM, JM, LM, 4) :: tausw, taulw
         real, dimension(LM, 4) :: dumtaubeam
         real, dimension(LM) :: dumasycl
         real, dimension(0:LM) :: dumtcldlyr
         real, dimension(0:LM) :: dumenn
         real :: taucir
         real, dimension(IM, JM, LM, 10) :: TAUDIAG
         real :: sunlit(IM * JM)
         real, dimension(IM, JM, LM) :: DELP, TAUSWICE, TAUSWLIQ, EMISS
         real, dimension(IM * JM, LM) :: EMISS2D, EMISSCOSP
         real :: isccp_emsfc_lw
         integer :: isccp_top_height
         integer :: isccp_top_height_direction

         ! icarus/isccp outputs
         real, dimension(IM * JM) :: ISCCP_totalcldarea
         real, dimension(IM * JM) :: ISCCP_meanptop
         real, dimension(IM * JM, 7, 7) :: fq_isccp
         real, dimension(IM, JM, 7, 7) :: fq_isccp3D
         real, dimension(IM * JM) :: ISCCP_meanalbedocld
         real, dimension(IM * JM) :: ISCCP_meanallskybrighttemp

         ! variables for lidar simulator
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         real, dimension(IM * JM, LM) :: lidar_pmol
         real, dimension(:, :, :), allocatable :: lidar_beta_tot, lidar_tau_tot
         real, dimension(IM * JM, PARASOL_NREFL) :: REFL

         ! variables for lidar and radar simulators
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         real, dimension(IM * JM, LM) :: ZLOCOSP
         real, dimension(IM * JM, LM) :: ZLO2D
         real :: K2
         integer :: use_mie_tables ! use a precomputed lookup table? yes=1,no=0,2=use first
         !  column everywhere
         integer :: use_gas_abs ! include gaseous absorption? yes=1,no=0
         integer :: do_ray ! calculate/output Rayleigh refl=1, not=0
         integer :: Nprmts_max_hydro ! Max number of parameters for hydrometeor size distributions
         integer :: Nprmts_max_aero ! Max number of parameters for aerosol size distributions
         real :: radar_freq
         integer :: lidar_ice_type !Ice particle shape in lidar calculations
         ! (0=ice-spheres;1=ice-non-spherical)

         ! radar and lidar simulator output -- the admittedly confusing lidar_ and radar_
         ! suffixes indicate the cosp return structures that the variables are copied out of
         real, dimension(:, :, :), allocatable :: radar_ze_tot
         real, dimension(:, :, :), allocatable :: radar_ze_tot_tmp
         logical, dimension(:, :, :), allocatable :: radar_ze_tot_mask
         integer, dimension(IM * JM, LM) :: ncvalid
         real, dimension(IM * JM, LM) :: radar_ze_tot_max
         real, dimension(IM * JM, LM) :: radar_ze_tot_mean

         real, dimension(IM * JM, dBZe_bins, LM) :: radar_cfad_ze
         real, dimension(IM * JM, LM) :: radar_lidar_only_freq_cloud
         real, dimension(IM * JM) :: radar_lidar_tcc
         real, dimension(IM * JM) :: land_mask

         real, dimension(IM * JM, SR_BINS, LM) :: lidar_cfad_sr ! CFAD of scattering ratio
         real, dimension(IM * JM, LM) :: lidar_lidarcld ! 3D "lidar" cloud fraction
         real, dimension(IM * JM, LIDAR_NCAT) :: lidar_cldlayer ! low, mid, high-level lidar cloud cover
         real, dimension(IM * JM, PARASOL_NREFL) :: lidar_parasolrefl ! mean parasol reflectance

         real, dimension(IM * JM, 7, MISR_N_CTH) :: fq_MISR
         real, dimension(IM * JM, MISR_N_CTH) :: MISR_dist_model_layertops
         real, dimension(IM * JM) :: MISR_meanztop
         real, dimension(IM * JM) :: MISR_cldarea

         real, dimension(IM * JM) :: MODIS_Cloud_Fraction_Total_Mean
         real, dimension(IM * JM) :: MODIS_Cloud_Fraction_Water_Mean
         real, dimension(IM * JM) :: MODIS_Cloud_Fraction_Ice_Mean
         real, dimension(IM * JM) :: MODIS_Cloud_Fraction_High_Mean
         real, dimension(IM * JM) :: MODIS_Cloud_Fraction_Mid_Mean
         real, dimension(IM * JM) :: MODIS_Cloud_Fraction_Low_Mean
         real, dimension(IM * JM) :: MODIS_Optical_Thickness_Total_Mean
         real, dimension(IM * JM) :: MODIS_Optical_Thickness_Water_Mean
         real, dimension(IM * JM) :: MODIS_Optical_Thickness_Ice_Mean
         real, dimension(IM * JM) :: MODIS_Optical_Thickness_Total_LogMean
         real, dimension(IM * JM) :: MODIS_Optical_Thickness_Water_LogMean
         real, dimension(IM * JM) :: MODIS_Optical_Thickness_Ice_LogMean
         real, dimension(IM * JM) :: MODIS_Cloud_Particle_Size_Water_Mean
         real, dimension(IM * JM) :: MODIS_Cloud_Particle_Size_Ice_Mean
         real, dimension(IM * JM) :: MODIS_Cloud_Top_Pressure_Total_Mean
         real, dimension(IM * JM) :: MODIS_Liquid_Water_Path_Mean
         real, dimension(IM * JM) :: MODIS_Ice_Water_Path_Mean
         real, dimension(IM * JM, 7, 7) :: MODIS_Optical_Thickness_vs_Cloud_Top_Pressure

         integer :: BEGSEG, ENDSEG

         type(SatSim_State), pointer :: self => null()
         type(SatSim_Wrap) :: wrap
         integer :: mapl_dims, vindex
         integer, pointer :: ungridded_dims(:)
         type(MAPL_VarSpec), pointer :: ExportSpec(:)
         real, pointer, dimension(:, :) :: ptr2d, ptr2d_new
         real, pointer, dimension(:, :) :: ptr_mask
         real, pointer, dimension(:, :, :) :: ptr3d, ptr3d_new
         type(ESMF_FieldBundle) :: bundle
         type(MAPL_MetaComp), pointer :: MAPL

         integer :: use_satsim, use_satsim_isccp, use_satsim_modis, use_satsim_lidar, use_satsim_radar, use_satsim_misr
         integer :: ncolumns

         character(len=ESMF_MAXSTR) :: GRIDNAME
         character(len=5) :: imchar
         character(len=2) :: dateline
         integer :: imsize, nn

         DEBUG_GC = .false.

         ! Get my MAPL_Generic state
         !--------------------------

         call MAPL_GetObjectFromGC(GC, MAPL, _RC)

         call MAPL_GetResource(MAPL, use_satsim, LABEL="USE_SATSIM:", default=0, _RC)

         call MAPL_GetResource(MAPL, use_satsim_isccp, LABEL="USE_SATSIM_ISCCP:", default=0, _RC)

         call MAPL_GetResource(MAPL, use_satsim_modis, LABEL="USE_SATSIM_MODIS:", default=0, _RC)

         call MAPL_GetResource(MAPL, use_satsim_radar, LABEL="USE_SATSIM_RADAR:", default=0, _RC)

         call MAPL_GetResource(MAPL, use_satsim_lidar, LABEL="USE_SATSIM_LIDAR:", default=0, _RC)

         call MAPL_GetResource(MAPL, use_satsim_misr, LABEL="USE_SATSIM_MISR:", default=0, _RC)

         call MAPL_GetResource(MAPL, GRIDNAME, 'AGCM_GRIDNAME:', _RC)
         GRIDNAME = AdjustL(GRIDNAME)
         nn = len_trim(GRIDNAME)
         dateline = GRIDNAME(nn - 1:nn)
         imchar = GRIDNAME(3:index(GRIDNAME, 'x') - 1)
         read(imchar, *) imsize
         if (dateline == 'CF') imsize = imsize * 4
         associate(default_Ncolumns => MIN(30, MAX(1, INT(4 * 5760 * 4 / imsize))))
            call MAPL_GetResource(MAPL, ncolumns, LABEL="SATSIM_NCOLUMNS:", default=default_Ncolumns, _RC)
         end associate

         call MAPL_GetResource(MAPL, Npoints_it, LABEL="SATSIM_POINTS_PER_ITERATION:", default=-999, _RC)

         allocate(frac_out(IM * JM, ncolumns, LM), _STAT)
         !allocate(      frac_outinv(IM*JM,NCOLUMNS,LM), _STAT)
         allocate(lidar_beta_tot(IM * JM, ncolumns, LM), _STAT)
         allocate(lidar_tau_tot(IM * JM, ncolumns, LM), _STAT)
         allocate(radar_ze_tot(IM * JM, ncolumns, LM), _STAT)
         allocate(radar_ze_tot_tmp(IM * JM, ncolumns, LM), _STAT)
         allocate(radar_ze_tot_mask(IM * JM, ncolumns, LM), _STAT)

         ! Pointers to imports
         !--------------------
         call MAPL_GetPointer(IMPORT, PLE, 'PLE', _RC)
         call MAPL_GetPointer(IMPORT, ZLE, 'ZLE', _RC)
         call MAPL_GetPointer(IMPORT, T, 'T', _RC)
         call MAPL_GetPointer(IMPORT, QV, 'QV', _RC)
         call MAPL_GetPointer(IMPORT, RDFL, 'RL', _RC)
         call MAPL_GetPointer(IMPORT, RDFI, 'RI', _RC)
         call MAPL_GetPointer(IMPORT, RDFR, 'RR', _RC)
         call MAPL_GetPointer(IMPORT, RDFS, 'RS', _RC)
         call MAPL_GetPointer(IMPORT, RDFG, 'RG', _RC)
         call MAPL_GetPointer(IMPORT, FCLD, 'FCLD', _RC)
         call MAPL_GetPointer(IMPORT, QLTOT, 'QL', _RC)
         call MAPL_GetPointer(IMPORT, QITOT, 'QI', _RC)
         call MAPL_GetPointer(IMPORT, QRTOT, 'QR', _RC)
         call MAPL_GetPointer(IMPORT, QSTOT, 'QS', _RC)
         call MAPL_GetPointer(IMPORT, QGTOT, 'QG', _RC)
         call MAPL_GetPointer(IMPORT, MCOSZ, 'MCOSZ', _RC)
         call MAPL_GetPointer(IMPORT, FRLAND, 'FRLAND', _RC)
         call MAPL_GetPointer(IMPORT, FROCEAN, 'FROCEAN', _RC)
         call MAPL_GetPointer(IMPORT, TS, 'TS', _RC)

         ! Pointers to Exports
         !--------------------

         ! ISCCP/Icarus

         call MAPL_GetPointer(EXPORT, ISCCP_CU_OA, 'ISCCP_CU_OA', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_CU_UA, 'ISCCP_CU_UA', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_CU_OB, 'ISCCP_CU_OB', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_CU_UB, 'ISCCP_CU_UB', _RC)

         call MAPL_GetPointer(EXPORT, ISCCP_STCU_OA, 'ISCCP_STCU_OA', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_STCU_UA, 'ISCCP_STCU_UA', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_STCU_OB, 'ISCCP_STCU_OB', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_STCU_UB, 'ISCCP_STCU_UB', _RC)

         call MAPL_GetPointer(EXPORT, ISCCP_ST_OA, 'ISCCP_ST_OA', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_ST_UA, 'ISCCP_ST_UA', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_ST_OB, 'ISCCP_ST_OB', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_ST_UB, 'ISCCP_ST_UB', _RC)

         call MAPL_GetPointer(EXPORT, ISCCP_ACU_OA, 'ISCCP_ACU_OA', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_ACU_UA, 'ISCCP_ACU_UA', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_ACU_OB, 'ISCCP_ACU_OB', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_ACU_UB, 'ISCCP_ACU_UB', _RC)

         call MAPL_GetPointer(EXPORT, ISCCP_AST_OA, 'ISCCP_AST_OA', _RC)

         call MAPL_GetPointer(EXPORT, ISCCP_AST_UA, 'ISCCP_AST_UA', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_AST_OB, 'ISCCP_AST_OB', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_AST_UB, 'ISCCP_AST_UB', _RC)

         call MAPL_GetPointer(EXPORT, ISCCP_NST_OA, 'ISCCP_NST_OA', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_NST_UA, 'ISCCP_NST_UA', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_NST_OB, 'ISCCP_NST_OB', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_NST_UB, 'ISCCP_NST_UB', _RC)

         call MAPL_GetPointer(EXPORT, ISCCP_CI_OA, 'ISCCP_CI_OA', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_CI_MA, 'ISCCP_CI_MA', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_CI_UA, 'ISCCP_CI_UA', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_CI_OB, 'ISCCP_CI_OB', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_CI_MB, 'ISCCP_CI_MB', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_CI_UB, 'ISCCP_CI_UB', _RC)

         call MAPL_GetPointer(EXPORT, ISCCP_CIST_OA, 'ISCCP_CIST_OA', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_CIST_MA, 'ISCCP_CIST_MA', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_CIST_UA, 'ISCCP_CIST_UA', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_CIST_OB, 'ISCCP_CIST_OB', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_CIST_MB, 'ISCCP_CIST_MB', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_CIST_UB, 'ISCCP_CIST_UB', _RC)

         call MAPL_GetPointer(EXPORT, ISCCP_CB_OA, 'ISCCP_CB_OA', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_CB_MA, 'ISCCP_CB_MA', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_CB_UA, 'ISCCP_CB_UA', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_CB_OB, 'ISCCP_CB_OB', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_CB_MB, 'ISCCP_CB_MB', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_CB_UB, 'ISCCP_CB_UB', _RC)

         call MAPL_GetPointer(EXPORT, ISCCP_SUBV1, 'ISCCP_SUBV1', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_SUBV2, 'ISCCP_SUBV2', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_SUBV3, 'ISCCP_SUBV3', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_SUBV4, 'ISCCP_SUBV4', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_SUBV5, 'ISCCP_SUBV5', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_SUBV6, 'ISCCP_SUBV6', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP_SUBV7, 'ISCCP_SUBV7', _RC)

         call MAPL_GetPointer(EXPORT, CLISCCP1, 'CLISCCP1', _RC)
         call MAPL_GetPointer(EXPORT, CLISCCP2, 'CLISCCP2', _RC)
         call MAPL_GetPointer(EXPORT, CLISCCP3, 'CLISCCP3', _RC)
         call MAPL_GetPointer(EXPORT, CLISCCP4, 'CLISCCP4', _RC)
         call MAPL_GetPointer(EXPORT, CLISCCP5, 'CLISCCP5', _RC)
         call MAPL_GetPointer(EXPORT, CLISCCP6, 'CLISCCP6', _RC)
         call MAPL_GetPointer(EXPORT, CLISCCP7, 'CLISCCP7', _RC)

         call MAPL_GetPointer(EXPORT, ISCCP1, 'ISCCP1', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP2, 'ISCCP2', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP3, 'ISCCP3', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP4, 'ISCCP4', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP5, 'ISCCP5', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP6, 'ISCCP6', _RC)
         call MAPL_GetPointer(EXPORT, ISCCP7, 'ISCCP7', _RC)

         call MAPL_GetPointer(EXPORT, TCLISCCP, 'TCLISCCP', _RC)
         call MAPL_GetPointer(EXPORT, CLTISCCP, 'CLTISCCP', _RC)
         call MAPL_GetPointer(EXPORT, PCTISCCP, 'PCTISCCP', _RC)
         call MAPL_GetPointer(EXPORT, CTPISCCP, 'CTPISCCP', _RC)
         call MAPL_GetPointer(EXPORT, ALBISCCP, 'ALBISCCP', _RC)
         call MAPL_GetPointer(EXPORT, TBISCCP, 'TBISCCP', _RC)
         call MAPL_GetPointer(EXPORT, SGFCLD, 'SGFCLD', _RC)

         ! LIDAR/CALIPSO

         call MAPL_GetPointer(EXPORT, RADARLTCC, 'RADARLTCC', _RC)
         call MAPL_GetPointer(EXPORT, CLCALIPSO2, 'CLCALIPSO2', _RC)
         call MAPL_GetPointer(EXPORT, LIDARPMOL, 'LIDARPMOL', _RC)
         call MAPL_GetPointer(EXPORT, LIDARPTOT, 'LIDARPTOT', _RC)
         call MAPL_GetPointer(EXPORT, LIDARTAUTOT, 'LIDARTAUTOT', _RC)
         call MAPL_GetPointer(EXPORT, RADARZETOT, 'RADARZETOT', _RC)
         call MAPL_GetPointer(EXPORT, CLCALIPSO, 'CLCALIPSO', _RC)
         call MAPL_GetPointer(EXPORT, CLLCALIPSO, 'CLLCALIPSO', _RC)
         call MAPL_GetPointer(EXPORT, CLMCALIPSO, 'CLMCALIPSO', _RC)
         call MAPL_GetPointer(EXPORT, CLHCALIPSO, 'CLHCALIPSO', _RC)
         call MAPL_GetPointer(EXPORT, CLTCALIPSO, 'CLTCALIPSO', _RC)
         call MAPL_GetPointer(EXPORT, PARASOLREFL0, 'PARASOLREFL0', _RC)
         call MAPL_GetPointer(EXPORT, PARASOLREFL1, 'PARASOLREFL1', _RC)
         call MAPL_GetPointer(EXPORT, PARASOLREFL2, 'PARASOLREFL2', _RC)
         call MAPL_GetPointer(EXPORT, PARASOLREFL3, 'PARASOLREFL3', _RC)
         call MAPL_GetPointer(EXPORT, PARASOLREFL4, 'PARASOLREFL4', _RC)
         call MAPL_GetPointer(EXPORT, PARASOLREFL5, 'PARASOLREFL5', _RC)

         call MAPL_GetPointer(EXPORT, CFADLIDARSR532_01, 'CFADLIDARSR532_01', _RC)
         call MAPL_GetPointer(EXPORT, CFADLIDARSR532_02, 'CFADLIDARSR532_02', _RC)
         call MAPL_GetPointer(EXPORT, CFADLIDARSR532_03, 'CFADLIDARSR532_03', _RC)
         call MAPL_GetPointer(EXPORT, CFADLIDARSR532_04, 'CFADLIDARSR532_04', _RC)
         call MAPL_GetPointer(EXPORT, CFADLIDARSR532_05, 'CFADLIDARSR532_05', _RC)
         call MAPL_GetPointer(EXPORT, CFADLIDARSR532_06, 'CFADLIDARSR532_06', _RC)
         call MAPL_GetPointer(EXPORT, CFADLIDARSR532_07, 'CFADLIDARSR532_07', _RC)
         call MAPL_GetPointer(EXPORT, CFADLIDARSR532_08, 'CFADLIDARSR532_08', _RC)
         call MAPL_GetPointer(EXPORT, CFADLIDARSR532_09, 'CFADLIDARSR532_09', _RC)
         call MAPL_GetPointer(EXPORT, CFADLIDARSR532_10, 'CFADLIDARSR532_10', _RC)
         call MAPL_GetPointer(EXPORT, CFADLIDARSR532_11, 'CFADLIDARSR532_11', _RC)
         call MAPL_GetPointer(EXPORT, CFADLIDARSR532_12, 'CFADLIDARSR532_12', _RC)
         call MAPL_GetPointer(EXPORT, CFADLIDARSR532_13, 'CFADLIDARSR532_13', _RC)
         call MAPL_GetPointer(EXPORT, CFADLIDARSR532_14, 'CFADLIDARSR532_14', _RC)
         call MAPL_GetPointer(EXPORT, CFADLIDARSR532_15, 'CFADLIDARSR532_15', _RC)

         ! RADAR/Cloudsat

         call MAPL_GetPointer(EXPORT, CLOUDSATCFAD01, 'CLOUDSATCFAD01', _RC)
         call MAPL_GetPointer(EXPORT, CLOUDSATCFAD02, 'CLOUDSATCFAD02', _RC)
         call MAPL_GetPointer(EXPORT, CLOUDSATCFAD03, 'CLOUDSATCFAD03', _RC)
         call MAPL_GetPointer(EXPORT, CLOUDSATCFAD04, 'CLOUDSATCFAD04', _RC)
         call MAPL_GetPointer(EXPORT, CLOUDSATCFAD05, 'CLOUDSATCFAD05', _RC)
         call MAPL_GetPointer(EXPORT, CLOUDSATCFAD06, 'CLOUDSATCFAD06', _RC)
         call MAPL_GetPointer(EXPORT, CLOUDSATCFAD07, 'CLOUDSATCFAD07', _RC)
         call MAPL_GetPointer(EXPORT, CLOUDSATCFAD08, 'CLOUDSATCFAD08', _RC)
         call MAPL_GetPointer(EXPORT, CLOUDSATCFAD09, 'CLOUDSATCFAD09', _RC)
         call MAPL_GetPointer(EXPORT, CLOUDSATCFAD10, 'CLOUDSATCFAD10', _RC)
         call MAPL_GetPointer(EXPORT, CLOUDSATCFAD11, 'CLOUDSATCFAD11', _RC)
         call MAPL_GetPointer(EXPORT, CLOUDSATCFAD12, 'CLOUDSATCFAD12', _RC)
         call MAPL_GetPointer(EXPORT, CLOUDSATCFAD13, 'CLOUDSATCFAD13', _RC)
         call MAPL_GetPointer(EXPORT, CLOUDSATCFAD14, 'CLOUDSATCFAD14', _RC)
         call MAPL_GetPointer(EXPORT, CLOUDSATCFAD15, 'CLOUDSATCFAD15', _RC)

         ! MODIS

         call MAPL_GetPointer(EXPORT, MDSCLDFRCTTL, 'MDSCLDFRCTTL', _RC)
         call MAPL_GetPointer(EXPORT, MDSCLDFRCWTR, 'MDSCLDFRCWTR', _RC)
         call MAPL_GetPointer(EXPORT, MDSCLDFRCH2O, 'MDSCLDFRCH2O', _RC)
         call MAPL_GetPointer(EXPORT, MDSCLDFRCICE, 'MDSCLDFRCICE', _RC)
         call MAPL_GetPointer(EXPORT, MDSCLDFRCHI, 'MDSCLDFRCHI', _RC)
         call MAPL_GetPointer(EXPORT, MDSCLDFRCMID, 'MDSCLDFRCMID', _RC)
         call MAPL_GetPointer(EXPORT, MDSCLDFRCLO, 'MDSCLDFRCLO', _RC)
         call MAPL_GetPointer(EXPORT, MDSOPTHCKTTL, 'MDSOPTHCKTTL', _RC)
         call MAPL_GetPointer(EXPORT, MDSOPTHCKWTR, 'MDSOPTHCKWTR', _RC)
         call MAPL_GetPointer(EXPORT, MDSOPTHCKH2O, 'MDSOPTHCKH2O', _RC)
         call MAPL_GetPointer(EXPORT, MDSOPTHCKICE, 'MDSOPTHCKICE', _RC)
         call MAPL_GetPointer(EXPORT, MDSOPTHCKTTLLG, 'MDSOPTHCKTTLLG', _RC)
         call MAPL_GetPointer(EXPORT, MDSOPTHCKWTRLG, 'MDSOPTHCKWTRLG', _RC)
         call MAPL_GetPointer(EXPORT, MDSOPTHCKH2OLG, 'MDSOPTHCKH2OLG', _RC)
         call MAPL_GetPointer(EXPORT, MDSOPTHCKICELG, 'MDSOPTHCKICELG', _RC)
         call MAPL_GetPointer(EXPORT, MDSCLDSZWTR, 'MDSCLDSZWTR', _RC)
         call MAPL_GetPointer(EXPORT, MDSCLDSZH20, 'MDSCLDSZH20', _RC)
         call MAPL_GetPointer(EXPORT, MDSCLDSZICE, 'MDSCLDSZICE', _RC)
         call MAPL_GetPointer(EXPORT, MDSCLDTOPPS, 'MDSCLDTOPPS', _RC)
         call MAPL_GetPointer(EXPORT, MDSWTRPATH, 'MDSWTRPATH', _RC)
         call MAPL_GetPointer(EXPORT, MDSH2OPATH, 'MDSH2OPATH', _RC)
         call MAPL_GetPointer(EXPORT, MDSICEPATH, 'MDSICEPATH', _RC)

         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST11, 'MDSTAUPRSHIST11', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST12, 'MDSTAUPRSHIST12', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST13, 'MDSTAUPRSHIST13', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST14, 'MDSTAUPRSHIST14', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST15, 'MDSTAUPRSHIST15', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST16, 'MDSTAUPRSHIST16', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST17, 'MDSTAUPRSHIST17', _RC)

         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST21, 'MDSTAUPRSHIST21', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST22, 'MDSTAUPRSHIST22', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST23, 'MDSTAUPRSHIST23', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST24, 'MDSTAUPRSHIST24', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST25, 'MDSTAUPRSHIST25', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST26, 'MDSTAUPRSHIST26', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST27, 'MDSTAUPRSHIST27', _RC)

         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST31, 'MDSTAUPRSHIST31', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST32, 'MDSTAUPRSHIST32', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST33, 'MDSTAUPRSHIST33', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST34, 'MDSTAUPRSHIST34', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST35, 'MDSTAUPRSHIST35', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST36, 'MDSTAUPRSHIST36', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST37, 'MDSTAUPRSHIST37', _RC)

         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST41, 'MDSTAUPRSHIST41', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST42, 'MDSTAUPRSHIST42', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST43, 'MDSTAUPRSHIST43', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST44, 'MDSTAUPRSHIST44', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST45, 'MDSTAUPRSHIST45', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST46, 'MDSTAUPRSHIST46', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST47, 'MDSTAUPRSHIST47', _RC)

         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST51, 'MDSTAUPRSHIST51', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST52, 'MDSTAUPRSHIST52', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST53, 'MDSTAUPRSHIST53', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST54, 'MDSTAUPRSHIST54', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST55, 'MDSTAUPRSHIST55', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST56, 'MDSTAUPRSHIST56', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST57, 'MDSTAUPRSHIST57', _RC)

         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST61, 'MDSTAUPRSHIST61', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST62, 'MDSTAUPRSHIST62', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST63, 'MDSTAUPRSHIST63', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST64, 'MDSTAUPRSHIST64', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST65, 'MDSTAUPRSHIST65', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST66, 'MDSTAUPRSHIST66', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST67, 'MDSTAUPRSHIST67', _RC)

         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST71, 'MDSTAUPRSHIST71', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST72, 'MDSTAUPRSHIST72', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST73, 'MDSTAUPRSHIST73', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST74, 'MDSTAUPRSHIST74', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST75, 'MDSTAUPRSHIST75', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST76, 'MDSTAUPRSHIST76', _RC)
         call MAPL_GetPointer(EXPORT, MDSTAUPRSHIST77, 'MDSTAUPRSHIST77', _RC)

         call MAPL_GetPointer(EXPORT, MISRMNCLDTP, 'MISRMNCLDTP', _RC)
         call MAPL_GetPointer(EXPORT, MISRCLDAREA, 'MISRCLDAREA', _RC)

         call MAPL_GetPointer(EXPORT, MISRLYRTP0, 'MISRLYRTP0', _RC)
         call MAPL_GetPointer(EXPORT, MISRLYRTP250, 'MISRLYRTP250', _RC)
         call MAPL_GetPointer(EXPORT, MISRLYRTP750, 'MISRLYRTP750', _RC)
         call MAPL_GetPointer(EXPORT, MISRLYRTP1250, 'MISRLYRTP1250', _RC)
         call MAPL_GetPointer(EXPORT, MISRLYRTP1750, 'MISRLYRTP1750', _RC)
         call MAPL_GetPointer(EXPORT, MISRLYRTP2250, 'MISRLYRTP2250', _RC)
         call MAPL_GetPointer(EXPORT, MISRLYRTP2750, 'MISRLYRTP2750', _RC)
         call MAPL_GetPointer(EXPORT, MISRLYRTP3500, 'MISRLYRTP3500', _RC)
         call MAPL_GetPointer(EXPORT, MISRLYRTP4500, 'MISRLYRTP4500', _RC)
         call MAPL_GetPointer(EXPORT, MISRLYRTP6000, 'MISRLYRTP6000', _RC)
         call MAPL_GetPointer(EXPORT, MISRLYRTP8000, 'MISRLYRTP8000', _RC)
         call MAPL_GetPointer(EXPORT, MISRLYRTP10000, 'MISRLYRTP10000', _RC)
         call MAPL_GetPointer(EXPORT, MISRLYRTP12000, 'MISRLYRTP12000', _RC)
         call MAPL_GetPointer(EXPORT, MISRLYRTP14000, 'MISRLYRTP14000', _RC)
         call MAPL_GetPointer(EXPORT, MISRLYRTP16000, 'MISRLYRTP16000', _RC)
         call MAPL_GetPointer(EXPORT, MISRLYRTP18000, 'MISRLYRTP18000', _RC)
         call MAPL_GetPointer(EXPORT, MISRFQ0, 'MISRFQ0', _RC)
         call MAPL_GetPointer(EXPORT, MISRFQ250, 'MISRFQ250', _RC)
         call MAPL_GetPointer(EXPORT, MISRFQ750, 'MISRFQ750', _RC)
         call MAPL_GetPointer(EXPORT, MISRFQ1250, 'MISRFQ1250', _RC)
         call MAPL_GetPointer(EXPORT, MISRFQ1750, 'MISRFQ1750', _RC)
         call MAPL_GetPointer(EXPORT, MISRFQ2250, 'MISRFQ2250', _RC)
         call MAPL_GetPointer(EXPORT, MISRFQ2750, 'MISRFQ2750', _RC)
         call MAPL_GetPointer(EXPORT, MISRFQ3500, 'MISRFQ3500', _RC)
         call MAPL_GetPointer(EXPORT, MISRFQ4500, 'MISRFQ4500', _RC)
         call MAPL_GetPointer(EXPORT, MISRFQ6000, 'MISRFQ6000', _RC)
         call MAPL_GetPointer(EXPORT, MISRFQ8000, 'MISRFQ8000', _RC)
         call MAPL_GetPointer(EXPORT, MISRFQ10000, 'MISRFQ10000', _RC)
         call MAPL_GetPointer(EXPORT, MISRFQ12000, 'MISRFQ12000', _RC)
         call MAPL_GetPointer(EXPORT, MISRFQ14000, 'MISRFQ14000', _RC)
         call MAPL_GetPointer(EXPORT, MISRFQ16000, 'MISRFQ16000', _RC)
         call MAPL_GetPointer(EXPORT, MISRFQ18000, 'MISRFQ18000', _RC)

         !  The follwing assignments are for variables needed for COSP and isccp/icarus
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         isccp_top_height = 1 !  1 = adjust top height using both a computed
         isccp_top_height_direction = 2 ! direction for finding atmosphere pressure level
         isccp_emsfc_lw = 0.999
         sunlit = 0.0
         icb = 58
         ict = 48

         ! anything that needs to be calculated in 3D is done so, to
         ! be converted to 2D and inverted later
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

         PLO = 0.5 * (PLE(:, :, 1:LM) + PLE(:, :, 0:LM - 1))

         ! The following are convective terms, which we don't want now
         CCA = 0.0

         ! the following is in lieu of
         !      where(FCLD.gt.0.)
         ! to avoid stupid high cwc, which makes for stupid high dtau_s, which
         ! crashes icarus
         where (FCLD > 0.01)
            CWC(:, :, :, 1) = MAX(QITOT / FCLD, 1.0e-12)
            CWC(:, :, :, 2) = MAX(QLTOT / FCLD, 1.0e-12)
            CWC(:, :, :, 3) = MAX(QRTOT / FCLD, 1.0e-12)
            CWC(:, :, :, 4) = MAX(QSTOT / FCLD, 1.0e-12)
            elsewhere
            CWC(:, :, :, 1) = 0.
            CWC(:, :, :, 2) = 0.
            CWC(:, :, :, 3) = 0.
            CWC(:, :, :, 4) = 0.
         end where

         ! delp in Pascals, reff in microns

#ifdef USE_MAPL_UNDEF
         where (RDFI == MAPL_UNDEF)
            RDFI = 36.e-6
         end where
         where (RDFL == MAPL_UNDEF)
            RDFL = 14.e-6
         end where
         where (RDFR == MAPL_UNDEF)
            RDFR = 50.e-6
         end where
         where (RDFS == MAPL_UNDEF)
            RDFS = 50.e-6
         end where
         where (RDFG == MAPL_UNDEF)
            RDFG = 50.e-6
         end where
#endif

         REFF(:, :, :, 1) = RDFI * 1.e6
         REFF(:, :, :, 2) = RDFL * 1.e6
         REFF(:, :, :, 3) = RDFR * 1.e6
         REFF(:, :, :, 4) = RDFS * 1.e6

         DELP = PLE(:, :, 1:LM) - PLE(:, :, 0:LM - 1)

         do i = 1, IM
            do j = 1, JM

               ! NOTE: mcosz is being passed in but not used since we aren't doing cloud scaling
               !       as the 0's at the end of inputs determine. If scaling is needed, mcosz might
               !       be the wrong cosine to use (i.e., not the mean). Also, dummies will need to
               !       be replaced if needed.
               call getvistau(LM, MCOSZ(i, j), DELP(i, j, :), FCLD(i, j, :), REFF(i, j, :, :), CWC(i, j, :, :), 0, 0,&
                    dumtaubeam(:, :), tausw(i, j, :, :), dumasycl(:))

               ! IR Taus: Here we calculate band 4 only. Again, dummies for unneeded arrays.
               call getirtau(4, LM, DELP(i, j, :), FCLD(i, j, :), REFF(i, j, :, :), CWC(i, j, :, :),&
                    taulw(i, j, :, :), dumtcldlyr(:), dumenn(:))

               do k = 1, LM

                  ! NOTE: If you wish to include both falling rain and falling snow to
                  !       DTAU_S and EMISS, use the lines below with all four
                  !       tausw and taulw returns

                  ! VIS
                  ! ---

                  !DTAU_S(i,j,k) = tausw(i,j,k,1)+tausw(i,j,k,2)+tausw(i,j,k,3)+tausw(i,j,k,4)
                  DTAU_S(i, j, k) = tausw(i, j, k, 1) + tausw(i, j, k, 2)

                  ! IR
                  ! --

                  !taucir = taulw(i,j,k,1) + taulw(i,j,k,2) + taulw(i,j,k,3) + taulw(i,j,k,4)
                  taucir = taulw(i, j, k, 1) + taulw(i, j, k, 2)

                  EMISS(i, j, k) = 1 - EXP(-1.0 * taucir)

               end do

            end do
         end do

         !  The follwing assignments are for variables needed for all simulators
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         Npoints = IM * JM
         isccp_overlap = 3 !  overlap type: 1=max, 2=rand, 3=max/rand
         ! used by SCOPS, COSP and icarus

         !  The follwing assignments are for variables needed for COSP
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         Naero = 1 ! Number of aerosol species (Not used)
         if (Npoints_it == -999) Npoints_it = Npoints ! Number of gridpoints processed in one iteration (probably not used)
         use_precipitation_fluxes = .false.
         use_reff = .true.
         surface_radar = 0 ! surface=1, spaceborne=0
         melt_lay = 0 ! melting layer model off=0, on=1
         time = 0 ! Time since start of run [days] - something should be done about this

         ! Switches in cfg need to be toggled to make things run right

         ! remember these are integers so + is logical or, not and

         if (use_satsim + use_satsim_isccp + use_satsim_modis > 0) then ! isccp needed for modis
            cfg%Lisccp_sim = .true.
         else
            cfg%Lisccp_sim = .false.
         end if

         if (use_satsim + use_satsim_modis > 0) then
            cfg%Lmodis_sim = .true.
         else
            cfg%Lmodis_sim = .false.
         end if

         if (use_satsim + use_satsim_radar > 0) then
            cfg%Lradar_sim = .true.
         else
            cfg%Lradar_sim = .false.
         end if

         if (use_satsim + use_satsim_lidar > 0) then
            cfg%Llidar_sim = .true.
            cfg%LCFADLIDARSR532 = .true.
         else
            cfg%Llidar_sim = .false.
            cfg%LCFADLIDARSR532 = .false.
         end if

         if (use_satsim + use_satsim_misr > 0) then
            cfg%Lmisr_sim = .true.
         else
            cfg%Lmisr_sim = .false.
         end if

         if (use_satsim /= 0 .or. (use_satsim_lidar + use_satsim_radar == 2)) then
            cfg%Lstats = .true.
         else
            cfg%Lstats = .false.
         end if

         ! RRTOV paramters (not used)
         Plat = 0
         Sat = 0
         Inst = 0
         Nchan = 0
         ZenAng = 0.0
         Ichan = 0
         SurfEm = 0.0
         co2 = 0.0
         ch4 = 0.0
         n2o = 0.0
         co = 0.0

         !  The follwing assignments are for variables needed for COSP and lidar_simulator
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         lidar_ice_type = 0 ! Ice particle shape in lidar calculations (0=ice-spheres ; 1=ice-non-spherical)

         !  The follwing assignments are for variables needed for COSP and radar_simulator
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         Nprmts_max_aero = 1 ! Max number of parameters for aerosol size distributions (Not used)
         Nprmts_max_hydro = 12 ! Max number of parameters for hydrometeor size distributions
         radar_freq = 94.0 ! CloudSat radar frequency (GHz)
         K2 = -1. ! |K|^2, -1=use frequency dependent default
         do_ray = 0 ! calculate/output Rayleigh refl=1, not=0
         use_gas_abs = 1 ! include gaseous absorption? yes=1,no=0
         use_mie_tables = 0 ! use a precomputed lookup table? yes=1,no=0
         QST = GEOS_QSAT(T, PLO / 100.)
         RH = (QV / QST) * 100.0 ! rel. humidity (in percent)
         land_mask = 0.0 !  needed for stats -- should be changed [0 - Ocean, 1 - Land]

         RHCOSP = reshape(RH(:, :, LM:1:-1), (/ IM * JM, LM /))
         PLOCOSP = reshape(PLO(:, :, LM:1:-1), (/ IM * JM, LM /))
         QVCOSP = reshape(QV(:, :, LM:1:-1), (/ IM * JM, LM /))
         CCACOSP = reshape(CCA(:, :, LM:1:-1), (/ IM * JM, LM /))
         dtau_sCOSP = reshape(DTAU_S(:, :, LM:1:-1), (/ IM * JM, LM /))
         TCOSP = reshape(T(:, :, LM:1:-1), (/ IM * JM, LM /))
         EMISSCOSP = reshape(EMISS(:, :, LM:1:-1), (/ IM * JM, LM /))
         FCLDCOSP = reshape(FCLD(:, :, LM:1:-1), (/ IM * JM, LM /))
         RDFLCOSP = reshape(RDFL(:, :, LM:1:-1), (/ IM * JM, LM /))
         RDFICOSP = reshape(RDFI(:, :, LM:1:-1), (/ IM * JM, LM /))
         RDFRCOSP = reshape(RDFR(:, :, LM:1:-1), (/ IM * JM, LM /))
         RDFSCOSP = reshape(RDFS(:, :, LM:1:-1), (/ IM * JM, LM /))
         RDFGCOSP = reshape(RDFG(:, :, LM:1:-1), (/ IM * JM, LM /))
         QLTOTCOSP = reshape(QLTOT(:, :, LM:1:-1), (/ IM * JM, LM /))
         QITOTCOSP = reshape(QITOT(:, :, LM:1:-1), (/ IM * JM, LM /))
         QRTOTCOSP = reshape(QRTOT(:, :, LM:1:-1), (/ IM * JM, LM /))
         QSTOTCOSP = reshape(QSTOT(:, :, LM:1:-1), (/ IM * JM, LM /))
         QGTOTCOSP = reshape(QGTOT(:, :, LM:1:-1), (/ IM * JM, LM /))
         PLECOSP = reshape(PLE(:, :, LM:0:-1), (/ IM * JM, LM + 1 /))
         MCOSZCOSP = reshape(MCOSZ, (/ IM * JM /))
         FRLANDCOSP = reshape(FRLAND, (/ IM * JM /))
         FROCEANCOSP = reshape(FROCEAN, (/ IM * JM /))
         TSCOSP = reshape(TS, (/ IM * JM /))

         ZLE2D = reshape(ZLE, (/ IM * JM, LM + 1 /))
         do k = 0, LM
            ZLE2D(:, k) = ZLE2D(:, k) - ZLE2D(:, LM)
         end do
         ZLECOSP = ZLE2D(:, LM:0:-1)
         ZLO2D = 0.5 * (ZLE2D(:, 0:LM - 1) + ZLE2D(:, 1:LM))
         ZLOCOSP = ZLO2D(:, LM:1:-1)

         where (MCOSZCOSP > 0.0) sunlit = 1.0
         where (FROCEANCOSP < 0.5) land_mask = 1.0

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! zero out exports
         ISCCP_totalcldarea = 0.0
         ISCCP_meanptop = 0.0
         ISCCP_meanalbedocld = 0.0
         ISCCP_meanallskybrighttemp = 0.0
         frac_out = 0.0
         fq_isccp = 0.0
         lidar_lidarcld = 0.0
         lidar_cldlayer = 0.0
         lidar_parasolrefl = 0.0
         lidar_cfad_sr = 0.0
         radar_lidar_tcc = 0.0
         radar_lidar_only_freq_cloud = 0.0
         radar_cfad_ze = 0.0
         radar_ze_tot = 0.0

         lidar_pmol = 0.0
         REFL = 0.0

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         call construct_cosp_gridbox(time &
              , time_bands &
              , radar_freq &
              , surface_radar &
              , use_mie_tables &
              , use_gas_abs &
              , do_ray &
              , melt_lay &
              , K2 &
              , Npoints &
              , LM &
              , ncolumns &
              , N_HYDRO &
              , Nprmts_max_hydro &
              , Naero &
              , Nprmts_max_aero &
              , Npoints_it &
              , lidar_ice_type &
              , isccp_top_height &
              , isccp_top_height_direction &
              , isccp_overlap &
              , isccp_emsfc_lw &
              , use_precipitation_fluxes &
              , use_reff &
              , Plat &
              , Sat &
              , Inst &
              , Nchan &
              , ZenAng &
              , Ichan &
              , SurfEm &
              , co2 &
              , ch4 &
              , n2o &
              , co &
              , gbx)

         call construct_cosp_subgrid(Npoints, ncolumns, LM, sgx)
         if (cfg%Lisccp_sim) call construct_cosp_isccp(cfg, Npoints, ncolumns, LM, isccp)
         !call construct_cosp_sghydro(Npoints,Ncolumns,LM,N_hydro,sghydro)
         if (cfg%Lradar_sim) call construct_cosp_sgradar(cfg, Npoints, ncolumns, LM, N_HYDRO, sgradar)
         if (cfg%Llidar_sim) call construct_cosp_sglidar(cfg, Npoints, ncolumns, LM, N_HYDRO, PARASOL_NREFL, sglidar)
         if (cfg%Lmodis_sim) call construct_cosp_modis(cfg, Npoints, modis)
         if (cfg%Lmisr_sim) call construct_cosp_misr(cfg, Npoints, misr)

         !   vgrid is use for lidar and radar simulation stats, but only if the stats are on
         ! a different grid from the model.  So mostly dummy params get passed`
         call construct_cosp_vgrid(gbx, 0, .false., .false., vgrid)

         if (cfg%Lradar_sim) call construct_cosp_radarstats(cfg, Npoints, ncolumns, LM, N_HYDRO, stradar)
         if (cfg%Llidar_sim) call construct_cosp_lidarstats(cfg, Npoints, ncolumns, LM, N_HYDRO, PARASOL_NREFL, stlidar)

         BEGSEG = 1
         ENDSEG = Npoints

         gbx%T = TCOSP(BEGSEG:ENDSEG, :)
         gbx%SH = QVCOSP(BEGSEG:ENDSEG, :)
         gbx%Q = RHCOSP(BEGSEG:ENDSEG, :)
         gbx%CCA = CCACOSP(BEGSEG:ENDSEG, :)
         gbx%TCA = FCLDCOSP(BEGSEG:ENDSEG, :)
         ! in cosp_types, gbx%ph is supposed to be half pressure levels, with the bottom as surface
         !   pressure -- see psfc in cosp_test
         gbx%ph = PLECOSP(BEGSEG:ENDSEG, 1:LM)
         gbx%p = PLOCOSP(BEGSEG:ENDSEG, :)
         gbx%psfc = PLECOSP(BEGSEG:ENDSEG, 0) ! surface pressure
         gbx%DTAU_S = dtau_sCOSP(BEGSEG:ENDSEG, :)
         gbx%dtau_c = 0.
         gbx%dem_s = EMISSCOSP(BEGSEG:ENDSEG, :)
         gbx%dem_c = 0.
         gbx%sunlit = sunlit(BEGSEG:ENDSEG)
         gbx%land = land_mask(BEGSEG:ENDSEG)
         gbx%skt = TSCOSP(BEGSEG:ENDSEG)
         gbx%mr_hydro(:, :, I_LSCLIQ) = QLTOTCOSP(BEGSEG:ENDSEG, :)
         gbx%mr_hydro(:, :, I_LSCICE) = QITOTCOSP(BEGSEG:ENDSEG, :)
         gbx%mr_hydro(:, :, I_LSRAIN) = QRTOTCOSP(BEGSEG:ENDSEG, :)
         gbx%mr_hydro(:, :, I_LSSNOW) = QSTOTCOSP(BEGSEG:ENDSEG, :)
         gbx%mr_hydro(:, :, I_LSGRPL) = QGTOTCOSP(BEGSEG:ENDSEG, :)
         gbx%mr_hydro(:, :, I_CVCLIQ) = 0.0
         gbx%mr_hydro(:, :, I_CVCICE) = 0.0
         gbx%mr_hydro(:, :, I_CVRAIN) = 0.0
         gbx%mr_hydro(:, :, I_CVSNOW) = 0.0
         gbx%REFF(:, :, I_LSCLIQ) = RDFLCOSP(BEGSEG:ENDSEG, :)
         gbx%REFF(:, :, I_LSCICE) = RDFICOSP(BEGSEG:ENDSEG, :)
         gbx%REFF(:, :, I_LSRAIN) = RDFRCOSP(BEGSEG:ENDSEG, :)
         gbx%REFF(:, :, I_LSSNOW) = RDFSCOSP(BEGSEG:ENDSEG, :)
         gbx%REFF(:, :, I_LSGRPL) = RDFGCOSP(BEGSEG:ENDSEG, :)
         gbx%zlev = ZLOCOSP(BEGSEG:ENDSEG, :)
         gbx%zlev_half = ZLECOSP(BEGSEG:ENDSEG, 1:LM)

         where (gbx%mr_hydro < 0.0) gbx%mr_hydro = 0.0

         if (MAPL_AM_I_Root() .and. DEBUG_GC) then

            write(*, *) 'ncolumns: ', ncolumns
            write(*, *) 'lm: ', LM
            write(*, *) 'npoints: ', Npoints
            write(*, *) 'npoints_it: ', Npoints_it

            !write(*,*) 'FCLD: ', FCLD

            !write(*,*) 'gbx%dtau_s: ', gbx%dtau_s
            !write(*,*) 'gbx%dem_s: ', gbx%dem_s

            !write(*,*) 'gbx%tca: ', gbx%tca
            !write(*,*) 'sgx%frac_out: ', sgx%frac_out

         end if

         call MAPL_TimerOn(STATE, "-COSP")

         call COSP(isccp_overlap, ncolumns, cfg, vgrid, gbx, sgx, sgradar, sglidar, isccp, misr, modis, stradar, &
              stlidar)

         call MAPL_TimerOff(STATE, "-COSP")

         ! Dims of sgx%frac_out (Npoints,Ncolumns,Nlevels)
         if (associated(SGFCLD)) then
            SGFCLD = reshape( &
                 sum(sgx%frac_out(:, :, LM:1:-1), 2) &
                 , (/ IM, JM, LM /)) &
                 / (1. * ncolumns)
         end if

         ! more outputs are available in the structures -- the current ones are those
         ! needed for CFMIP

         if (use_satsim_lidar + use_satsim > 0) then

            lidar_parasolrefl(BEGSEG:ENDSEG, :) = stlidar%parasolrefl
            lidar_cldlayer(BEGSEG:ENDSEG, :) = stlidar%cldlayer
            lidar_lidarcld(BEGSEG:ENDSEG, :) = stlidar%lidarcld
            lidar_cfad_sr(BEGSEG:ENDSEG, :, :) = stlidar%cfad_sr

            ! used for diagnostics
            lidar_pmol(BEGSEG:ENDSEG, :) = sglidar%beta_mol
            lidar_beta_tot(BEGSEG:ENDSEG, :, :) = sglidar%beta_tot
            lidar_tau_tot(BEGSEG:ENDSEG, :, :) = sglidar%tau_tot

         end if

         if (use_satsim_radar + use_satsim > 0) then
            radar_ze_tot(BEGSEG:ENDSEG, :, :) = sgradar%ze_tot
            radar_lidar_only_freq_cloud(BEGSEG:ENDSEG, :) = stradar%lidar_only_freq_cloud
            radar_lidar_tcc(BEGSEG:ENDSEG) = stradar%radar_lidar_tcc
            radar_cfad_ze(BEGSEG:ENDSEG, :, :) = stradar%cfad_ze
         end if

         if (use_satsim + use_satsim_isccp > 0) then
            ! cosp reverses the pressure levls of the isccp matrix for some reason
            fq_isccp(BEGSEG:ENDSEG, :, :) = isccp%fq_isccp(:, :, 7:1:-1)
            ISCCP_totalcldarea(BEGSEG:ENDSEG) = isccp%totalcldarea
            ISCCP_meanptop(BEGSEG:ENDSEG) = isccp%meanptop
            ISCCP_meanalbedocld(BEGSEG:ENDSEG) = isccp%meanalbedocld
            ISCCP_meanallskybrighttemp(BEGSEG:ENDSEG) = isccp%meantb
            fq_isccp3D = reshape(fq_isccp, (/ IM, JM, 7, 7 /))

#ifdef USE_MAPL_UNDEF

            where (fq_isccp3D < -10.0) fq_isccp3D = MAPL_UNDEF
            where (fq_isccp3D > 10.0) fq_isccp3D = MAPL_UNDEF
            where (ISCCP_totalcldarea < 0.0) ISCCP_totalcldarea = MAPL_UNDEF
            where (ISCCP_totalcldarea > 1.0) ISCCP_totalcldarea = MAPL_UNDEF
            where (ISCCP_meanptop <= 0.0) ISCCP_meanptop = MAPL_UNDEF
            where (ISCCP_meanalbedocld < 0.0) ISCCP_meanalbedocld = MAPL_UNDEF
            where (ISCCP_meanalbedocld > 1.0) ISCCP_meanalbedocld = MAPL_UNDEF
            where (ISCCP_meanallskybrighttemp < 0.0) ISCCP_meanallskybrighttemp = MAPL_UNDEF

#endif

         end if ! ( USE_SATSIM + USE_SATSIM_ISCCP > 0 )

         if (use_satsim_misr + use_satsim > 0) then
            MISR_meanztop(BEGSEG:ENDSEG) = misr%MISR_meanztop
            fq_MISR(BEGSEG:ENDSEG, :, :) = misr%fq_MISR
            MISR_cldarea(BEGSEG:ENDSEG) = misr%MISR_cldarea
            MISR_dist_model_layertops(BEGSEG:ENDSEG, :) = misr%MISR_dist_model_layertops
         end if

         if (use_satsim_modis + use_satsim > 0) then
            MODIS_Cloud_Fraction_Total_Mean(BEGSEG:ENDSEG) = modis%Cloud_Fraction_Total_Mean
            MODIS_Cloud_Fraction_Water_Mean(BEGSEG:ENDSEG) = modis%Cloud_Fraction_Water_Mean
            MODIS_Cloud_Fraction_Ice_Mean(BEGSEG:ENDSEG) = modis%Cloud_Fraction_Ice_Mean
            MODIS_Cloud_Fraction_High_Mean(BEGSEG:ENDSEG) = modis%Cloud_Fraction_High_Mean
            MODIS_Cloud_Fraction_Mid_Mean(BEGSEG:ENDSEG) = modis%Cloud_Fraction_Mid_Mean
            MODIS_Cloud_Fraction_Low_Mean(BEGSEG:ENDSEG) = modis%Cloud_Fraction_Low_Mean
            MODIS_Optical_Thickness_Total_Mean(BEGSEG:ENDSEG) = modis%Optical_Thickness_Total_Mean
            MODIS_Optical_Thickness_Water_Mean(BEGSEG:ENDSEG) = modis%Optical_Thickness_Water_Mean
            MODIS_Optical_Thickness_Ice_Mean(BEGSEG:ENDSEG) = modis%Optical_Thickness_Ice_Mean
            MODIS_Optical_Thickness_Total_LogMean(BEGSEG:ENDSEG) = modis%Optical_Thickness_Total_LogMean
            MODIS_Optical_Thickness_Water_LogMean(BEGSEG:ENDSEG) = modis%Optical_Thickness_Water_LogMean
            MODIS_Optical_Thickness_Ice_LogMean(BEGSEG:ENDSEG) = modis%Optical_Thickness_Ice_LogMean
            MODIS_Cloud_Particle_Size_Water_Mean(BEGSEG:ENDSEG) = modis%Cloud_Particle_Size_Water_Mean
            MODIS_Cloud_Particle_Size_Ice_Mean(BEGSEG:ENDSEG) = modis%Cloud_Particle_Size_Ice_Mean
            MODIS_Cloud_Top_Pressure_Total_Mean(BEGSEG:ENDSEG) = modis%Cloud_Top_Pressure_Total_Mean
            MODIS_Liquid_Water_Path_Mean(BEGSEG:ENDSEG) = modis%Liquid_Water_Path_Mean
            MODIS_Ice_Water_Path_Mean(BEGSEG:ENDSEG) = modis%Ice_Water_Path_Mean
            MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(BEGSEG:ENDSEG, :, :) = modis%&
                 Optical_Thickness_vs_Cloud_Top_Pressure
         end if

#ifdef USE_MAPL_UNDEF

         if (use_satsim_modis + use_satsim > 0) then
            ! setting to undef since 0 implies no clouds
            where (MODIS_Optical_Thickness_Total_Mean <= 0.0) MODIS_Optical_Thickness_Total_Mean = MAPL_UNDEF
            where (MODIS_Optical_Thickness_Water_Mean <= 0.0) MODIS_Optical_Thickness_Water_Mean = MAPL_UNDEF
            where (MODIS_Optical_Thickness_Ice_Mean <= 0.0) MODIS_Optical_Thickness_Ice_Mean = MAPL_UNDEF
            where (MODIS_Optical_Thickness_Total_LogMean <= 0.0) MODIS_Optical_Thickness_Total_LogMean = MAPL_UNDEF
            where (MODIS_Optical_Thickness_Water_LogMean <= 0.0) MODIS_Optical_Thickness_Water_LogMean = MAPL_UNDEF
            where (MODIS_Optical_Thickness_Ice_LogMean <= 0.0) MODIS_Optical_Thickness_Ice_LogMean = MAPL_UNDEF
            where (MODIS_Cloud_Particle_Size_Water_Mean <= 0.0) MODIS_Cloud_Particle_Size_Water_Mean = MAPL_UNDEF
            where (MODIS_Cloud_Particle_Size_Ice_Mean <= 0.0) MODIS_Cloud_Particle_Size_Ice_Mean = MAPL_UNDEF

            where (MODIS_Cloud_Fraction_Total_Mean == R_UNDEF) MODIS_Cloud_Fraction_Total_Mean = MAPL_UNDEF
            where (MODIS_Cloud_Fraction_Water_Mean == R_UNDEF) MODIS_Cloud_Fraction_Water_Mean = MAPL_UNDEF
            where (MODIS_Cloud_Fraction_Ice_Mean == R_UNDEF) MODIS_Cloud_Fraction_Ice_Mean = MAPL_UNDEF
            where (MODIS_Cloud_Fraction_High_Mean == R_UNDEF) MODIS_Cloud_Fraction_High_Mean = MAPL_UNDEF
            where (MODIS_Cloud_Fraction_Mid_Mean == R_UNDEF) MODIS_Cloud_Fraction_Mid_Mean = MAPL_UNDEF
            where (MODIS_Cloud_Fraction_Low_Mean == R_UNDEF) MODIS_Cloud_Fraction_Low_Mean = MAPL_UNDEF
            where (MODIS_Optical_Thickness_Total_Mean == R_UNDEF) MODIS_Optical_Thickness_Total_Mean = MAPL_UNDEF
            where (MODIS_Optical_Thickness_Water_Mean == R_UNDEF) MODIS_Optical_Thickness_Water_Mean = MAPL_UNDEF
            where (MODIS_Optical_Thickness_Ice_Mean == R_UNDEF) MODIS_Optical_Thickness_Ice_Mean = MAPL_UNDEF
            where (MODIS_Optical_Thickness_Total_LogMean == R_UNDEF) MODIS_Optical_Thickness_Total_LogMean = MAPL_UNDEF
            where (MODIS_Optical_Thickness_Water_LogMean == R_UNDEF) MODIS_Optical_Thickness_Water_LogMean = MAPL_UNDEF
            where (MODIS_Optical_Thickness_Ice_LogMean == R_UNDEF) MODIS_Optical_Thickness_Ice_LogMean = MAPL_UNDEF
            where (MODIS_Cloud_Particle_Size_Water_Mean == R_UNDEF) MODIS_Cloud_Particle_Size_Water_Mean = MAPL_UNDEF
            where (MODIS_Cloud_Particle_Size_Ice_Mean == R_UNDEF) MODIS_Cloud_Particle_Size_Ice_Mean = MAPL_UNDEF
            where (MODIS_Cloud_Top_Pressure_Total_Mean == R_UNDEF) MODIS_Cloud_Top_Pressure_Total_Mean = MAPL_UNDEF
            where (MODIS_Cloud_Top_Pressure_Total_Mean <= 0.0) MODIS_Cloud_Top_Pressure_Total_Mean = MAPL_UNDEF
            where (MODIS_Liquid_Water_Path_Mean == R_UNDEF) MODIS_Liquid_Water_Path_Mean = MAPL_UNDEF
            where (MODIS_Ice_Water_Path_Mean == R_UNDEF) MODIS_Ice_Water_Path_Mean = MAPL_UNDEF
            where (MODIS_Optical_Thickness_vs_Cloud_Top_Pressure == R_UNDEF) &
                 MODIS_Optical_Thickness_vs_Cloud_Top_Pressure = MAPL_UNDEF
         end if

         if (use_satsim_lidar + use_satsim > 0) then
            where (lidar_cfad_sr == R_UNDEF) lidar_cfad_sr = MAPL_UNDEF
            where (lidar_lidarcld > 1.0) lidar_lidarcld = MAPL_UNDEF
            where (lidar_lidarcld < 0.0) lidar_lidarcld = MAPL_UNDEF
            where (lidar_cldlayer == R_UNDEF) lidar_cldlayer = MAPL_UNDEF
            where (lidar_parasolrefl == R_UNDEF) lidar_parasolrefl = MAPL_UNDEF
         end if

         if (use_satsim_radar + use_satsim > 0) then
            where (radar_ze_tot == R_UNDEF) radar_ze_tot = MAPL_UNDEF
            where (radar_lidar_only_freq_cloud == R_UNDEF) radar_lidar_only_freq_cloud = MAPL_UNDEF
            where (radar_lidar_tcc == R_UNDEF) radar_lidar_tcc = MAPL_UNDEF
            where (radar_cfad_ze == R_UNDEF) radar_cfad_ze = MAPL_UNDEF
         end if

         if (use_satsim_misr + use_satsim > 0) then
            where (MISR_meanztop == R_UNDEF) MISR_meanztop = MAPL_UNDEF
            where (MISR_meanztop <= 0.0) MISR_meanztop = MAPL_UNDEF
            where (fq_MISR == R_UNDEF) fq_MISR = MAPL_UNDEF
            where (MISR_cldarea == R_UNDEF) MISR_cldarea = MAPL_UNDEF
            where (MISR_dist_model_layertops == R_UNDEF) MISR_dist_model_layertops = MAPL_UNDEF
         end if
#endif

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! ISCCP exports
         !-------------------------------------------------
         !-------------------------------------------------

         ! 4 Cumulus (CU) subcategories
         !    fq_isccp3D(2:3,6:7)
         !-------------------------------------------------

         if (associated(ISCCP_CU_OA)) then
            ISCCP_CU_OA = fq_isccp3D(:, :, 2, 6)
         end if

         if (associated(ISCCP_CU_OB)) then
            ISCCP_CU_OB = fq_isccp3D(:, :, 3, 6)
         end if

         if (associated(ISCCP_CU_UA)) then
            ISCCP_CU_UA = fq_isccp3D(:, :, 2, 7)
         end if

         if (associated(ISCCP_CU_UB)) then
            ISCCP_CU_UB = fq_isccp3D(:, :, 3, 7)
         end if

         ! 4 Stratocumulus STCU subcategories
         !    ISCCP_isccp3D(4:5,6:7)
         !-------------------------------------------------

         if (associated(ISCCP_STCU_OA)) then
            ISCCP_STCU_OA = fq_isccp3D(:, :, 4, 6)
         end if

         if (associated(ISCCP_STCU_OB)) then
            ISCCP_STCU_OB = fq_isccp3D(:, :, 5, 6)
         end if

         if (associated(ISCCP_STCU_UA)) then
            ISCCP_STCU_UA = fq_isccp3D(:, :, 4, 7)
         end if

         if (associated(ISCCP_STCU_UB)) then
            ISCCP_STCU_UB = fq_isccp3D(:, :, 5, 7)
         end if

         ! 4 Stratus ST subcategories
         !    fq_isccp3D(6:7,6:7)
         !-------------------------------------------------

         if (associated(ISCCP_ST_OA)) then
            ISCCP_ST_OA = fq_isccp3D(:, :, 6, 6)
         end if

         if (associated(ISCCP_ST_OB)) then
            ISCCP_ST_OB = fq_isccp3D(:, :, 7, 6)
         end if

         if (associated(ISCCP_ST_UA)) then
            ISCCP_ST_UA = fq_isccp3D(:, :, 6, 7)
         end if

         if (associated(ISCCP_ST_UB)) then
            ISCCP_ST_UB = fq_isccp3D(:, :, 7, 7)
         end if

         ! 4 Altocumulus ACU subcategories
         !    fq_isccp3D(2:3,4:5)
         !-------------------------------------------------

         if (associated(ISCCP_ACU_OA)) then
            ISCCP_ACU_OA = fq_isccp3D(:, :, 2, 4)
         end if

         if (associated(ISCCP_ACU_OB)) then
            ISCCP_ACU_OB = fq_isccp3D(:, :, 3, 4)
         end if

         if (associated(ISCCP_ACU_UA)) then
            ISCCP_ACU_UA = fq_isccp3D(:, :, 2, 5)
         end if

         if (associated(ISCCP_ACU_UB)) then
            ISCCP_ACU_UB = fq_isccp3D(:, :, 3, 5)
         end if

         ! 4 Altostratus AST subcategories
         !    fq_isccp3D(4:5,4:5)
         !-------------------------------------------------

         if (associated(ISCCP_AST_OA)) then
            ISCCP_AST_OA = fq_isccp3D(:, :, 4, 4)
         end if

         if (associated(ISCCP_AST_OB)) then
            ISCCP_AST_OB = fq_isccp3D(:, :, 5, 4)
         end if

         if (associated(ISCCP_AST_UA)) then
            ISCCP_AST_UA = fq_isccp3D(:, :, 4, 5)
         end if

         if (associated(ISCCP_AST_UB)) then
            ISCCP_AST_UB = fq_isccp3D(:, :, 5, 5)
         end if

         ! 4 Nimbostratus NST subcategories
         !    fq_isccp3D(6:7,4:5)
         !-------------------------------------------------

         if (associated(ISCCP_NST_OA)) then
            ISCCP_NST_OA = fq_isccp3D(:, :, 6, 4)
         end if

         if (associated(ISCCP_NST_OB)) then
            ISCCP_NST_OB = fq_isccp3D(:, :, 7, 4)
         end if

         if (associated(ISCCP_NST_UA)) then
            ISCCP_NST_UA = fq_isccp3D(:, :, 6, 5)
         end if

         if (associated(ISCCP_NST_UB)) then
            ISCCP_NST_UB = fq_isccp3D(:, :, 7, 5)
         end if

         ! 6 Cirrus CI subcategories
         !    fq_isccp3D(2:3,1:3)
         !-------------------------------------------------

         if (associated(ISCCP_CI_OA)) then
            ISCCP_CI_OA = fq_isccp3D(:, :, 2, 1)
         end if

         if (associated(ISCCP_CI_OB)) then
            ISCCP_CI_OB = fq_isccp3D(:, :, 3, 1)
         end if

         if (associated(ISCCP_CI_MA)) then
            ISCCP_CI_MA = fq_isccp3D(:, :, 2, 2)
         end if

         if (associated(ISCCP_CI_MB)) then
            ISCCP_CI_MB = fq_isccp3D(:, :, 3, 2)
         end if

         if (associated(ISCCP_CI_UA)) then
            ISCCP_CI_UA = fq_isccp3D(:, :, 2, 3)
         end if

         if (associated(ISCCP_CI_UB)) then
            ISCCP_CI_UB = fq_isccp3D(:, :, 3, 3)
         end if

         ! 6 Cirrostratus CIST subcategories
         !    fq_isccp3D(4:5,1:3)
         !-------------------------------------------------

         if (associated(ISCCP_CIST_OA)) then
            ISCCP_CIST_OA = fq_isccp3D(:, :, 4, 1)
         end if

         if (associated(ISCCP_CIST_OB)) then
            ISCCP_CIST_OB = fq_isccp3D(:, :, 5, 1)
         end if

         if (associated(ISCCP_CIST_MA)) then
            ISCCP_CIST_MA = fq_isccp3D(:, :, 4, 2)
         end if

         if (associated(ISCCP_CIST_MB)) then
            ISCCP_CIST_MB = fq_isccp3D(:, :, 5, 2)
         end if

         if (associated(ISCCP_CIST_UA)) then
            ISCCP_CIST_UA = fq_isccp3D(:, :, 4, 3)
         end if

         if (associated(ISCCP_CIST_UB)) then
            ISCCP_CIST_UB = fq_isccp3D(:, :, 5, 3)
         end if

         ! 6 Cumulonimbus CB subcategories
         !    fq_isccp3D(6:7,1:3)
         !-------------------------------------------------

         if (associated(ISCCP_CB_OA)) then
            ISCCP_CB_OA = fq_isccp3D(:, :, 6, 1)
         end if

         if (associated(ISCCP_CB_OB)) then
            ISCCP_CB_OB = fq_isccp3D(:, :, 7, 1)
         end if

         if (associated(ISCCP_CB_MA)) then
            ISCCP_CB_MA = fq_isccp3D(:, :, 6, 2)
         end if

         if (associated(ISCCP_CB_MB)) then
            ISCCP_CB_MB = fq_isccp3D(:, :, 7, 2)
         end if

         if (associated(ISCCP_CB_UA)) then
            ISCCP_CB_UA = fq_isccp3D(:, :, 6, 3)
         end if

         if (associated(ISCCP_CB_UB)) then
            ISCCP_CB_UB = fq_isccp3D(:, :, 7, 3)
         end if

         ! 7 Subvisible/Subdetection (SUBV) subcategories
         !    ISCCP_isccp3D(1,1:7)
         !-------------------------------------------------

         if (associated(ISCCP_SUBV1)) then
            ISCCP_SUBV1 = fq_isccp3D(:, :, 1, 1)
         end if

         if (associated(ISCCP_SUBV2)) then
            ISCCP_SUBV2 = fq_isccp3D(:, :, 1, 2)
         end if

         if (associated(ISCCP_SUBV3)) then
            ISCCP_SUBV3 = fq_isccp3D(:, :, 1, 3)
         end if

         if (associated(ISCCP_SUBV4)) then
            ISCCP_SUBV4 = fq_isccp3D(:, :, 1, 4)
         end if

         if (associated(ISCCP_SUBV5)) then
            ISCCP_SUBV5 = fq_isccp3D(:, :, 1, 5)
         end if

         if (associated(ISCCP_SUBV6)) then
            ISCCP_SUBV6 = fq_isccp3D(:, :, 1, 6)
         end if

         if (associated(ISCCP_SUBV7)) then
            ISCCP_SUBV7 = fq_isccp3D(:, :, 1, 7)
         end if

         ! 7 pressure levels of all thicknesses
         !    ISCCP_isccp3D(:,:,:,1:7)
         !-------------------------------------------------

         if (associated(CLISCCP1)) then
            CLISCCP1 = fq_isccp3D(:, :, :, 1)
         end if

         if (associated(CLISCCP2)) then
            CLISCCP2 = fq_isccp3D(:, :, :, 2)
         end if

         if (associated(CLISCCP3)) then
            CLISCCP3 = fq_isccp3D(:, :, :, 3)
         end if

         if (associated(CLISCCP4)) then
            CLISCCP4 = fq_isccp3D(:, :, :, 4)
         end if

         if (associated(CLISCCP5)) then
            CLISCCP5 = fq_isccp3D(:, :, :, 5)
         end if

         if (associated(CLISCCP6)) then
            CLISCCP6 = fq_isccp3D(:, :, :, 6)
         end if

         if (associated(CLISCCP7)) then
            CLISCCP7 = fq_isccp3D(:, :, :, 7)
         end if

         if (associated(ISCCP1)) then
            ISCCP1 = fq_isccp3D(:, :, :, 1)
         end if

         if (associated(ISCCP2)) then
            ISCCP2 = fq_isccp3D(:, :, :, 2)
         end if

         if (associated(ISCCP3)) then
            ISCCP3 = fq_isccp3D(:, :, :, 3)
         end if

         if (associated(ISCCP4)) then
            ISCCP4 = fq_isccp3D(:, :, :, 4)
         end if

         if (associated(ISCCP5)) then
            ISCCP5 = fq_isccp3D(:, :, :, 5)
         end if

         if (associated(ISCCP6)) then
            ISCCP6 = fq_isccp3D(:, :, :, 6)
         end if

         if (associated(ISCCP7)) then
            ISCCP7 = fq_isccp3D(:, :, :, 7)
         end if

         if (associated(CLTISCCP)) then
            CLTISCCP = reshape(ISCCP_totalcldarea, (/ IM, JM /))
         end if

         if (associated(TCLISCCP)) then
            TCLISCCP = reshape(ISCCP_totalcldarea, (/ IM, JM /))
         end if

         if (associated(PCTISCCP)) then
            PCTISCCP = reshape(ISCCP_meanptop, (/ IM, JM /))
         end if

         if (associated(CTPISCCP)) then
            CTPISCCP = reshape(ISCCP_meanptop, (/ IM, JM /))
         end if

         if (associated(ALBISCCP)) then
            ALBISCCP = reshape(ISCCP_meanalbedocld, (/ IM, JM /))
         end if

         if (associated(TBISCCP)) then
            TBISCCP = reshape(ISCCP_meanallskybrighttemp, (/ IM, JM /))
         end if

         ! CALIPSO/lidar exports
         !-------------------------------------------------
         !-------------------------------------------------

         if (associated(CLCALIPSO)) then
            CLCALIPSO = reshape(lidar_lidarcld(:, LM:1:-1), (/IM, JM, LM /))
         end if

         if (associated(CLLCALIPSO)) then
            CLLCALIPSO = reshape(lidar_cldlayer(:, 1), (/IM, JM /)) ! low
         end if

         if (associated(CLMCALIPSO)) then
            CLMCALIPSO = reshape(lidar_cldlayer(:, 2), (/IM, JM /)) ! middle
         end if

         if (associated(CLHCALIPSO)) then
            CLHCALIPSO = reshape(lidar_cldlayer(:, 3), (/IM, JM /)) ! high
         end if

         if (associated(CLTCALIPSO)) then
            CLTCALIPSO = reshape(lidar_cldlayer(:, 4), (/IM, JM /)) ! total
         end if

         if (associated(PARASOLREFL0)) then
            PARASOLREFL0 = reshape(lidar_parasolrefl, (/IM, JM, PARASOL_NREFL /))
         end if

         if (associated(PARASOLREFL1)) then
            PARASOLREFL1 = reshape(lidar_parasolrefl(:, 1), (/IM, JM /))
         end if

         if (associated(PARASOLREFL2)) then
            PARASOLREFL2 = reshape(lidar_parasolrefl(:, 2), (/IM, JM /))
         end if

         if (associated(PARASOLREFL3)) then
            PARASOLREFL3 = reshape(lidar_parasolrefl(:, 3), (/IM, JM /))
         end if

         if (associated(PARASOLREFL4)) then
            PARASOLREFL4 = reshape(lidar_parasolrefl(:, 4), (/IM, JM /))
         end if

         if (associated(PARASOLREFL5)) then
            PARASOLREFL5 = reshape(lidar_parasolrefl(:, 5), (/IM, JM /))
         end if

         if (associated(CFADLIDARSR532_01)) then
            CFADLIDARSR532_01 = reshape(lidar_cfad_sr(:, 1, LM:1:-1), (/IM, JM, LM /))
         end if

         if (associated(CFADLIDARSR532_02)) then
            CFADLIDARSR532_02 = reshape(lidar_cfad_sr(:, 2, LM:1:-1), (/IM, JM, LM /))
         end if

         if (associated(CFADLIDARSR532_03)) then
            CFADLIDARSR532_03 = reshape(lidar_cfad_sr(:, 3, LM:1:-1), (/IM, JM, LM /))
         end if

         if (associated(CFADLIDARSR532_04)) then
            CFADLIDARSR532_04 = reshape(lidar_cfad_sr(:, 4, LM:1:-1), (/IM, JM, LM /))
         end if

         if (associated(CFADLIDARSR532_05)) then
            CFADLIDARSR532_05 = reshape(lidar_cfad_sr(:, 5, LM:1:-1), (/IM, JM, LM /))
         end if

         if (associated(CFADLIDARSR532_06)) then
            CFADLIDARSR532_06 = reshape(lidar_cfad_sr(:, 6, LM:1:-1), (/IM, JM, LM /))
         end if

         if (associated(CFADLIDARSR532_07)) then
            CFADLIDARSR532_07 = reshape(lidar_cfad_sr(:, 7, LM:1:-1), (/IM, JM, LM /))
         end if

         if (associated(CFADLIDARSR532_08)) then
            CFADLIDARSR532_08 = reshape(lidar_cfad_sr(:, 8, LM:1:-1), (/IM, JM, LM /))
         end if

         if (associated(CFADLIDARSR532_09)) then
            CFADLIDARSR532_09 = reshape(lidar_cfad_sr(:, 9, LM:1:-1), (/IM, JM, LM /))
         end if

         if (associated(CFADLIDARSR532_10)) then
            CFADLIDARSR532_10 = reshape(lidar_cfad_sr(:, 10, LM:1:-1), (/IM, JM, LM /))
         end if

         if (associated(CFADLIDARSR532_11)) then
            CFADLIDARSR532_11 = reshape(lidar_cfad_sr(:, 11, LM:1:-1), (/IM, JM, LM /))
         end if

         if (associated(CFADLIDARSR532_12)) then
            CFADLIDARSR532_12 = reshape(lidar_cfad_sr(:, 12, LM:1:-1), (/IM, JM, LM /))
         end if

         if (associated(CFADLIDARSR532_13)) then
            CFADLIDARSR532_13 = reshape(lidar_cfad_sr(:, 13, LM:1:-1), (/IM, JM, LM /))
         end if

         if (associated(CFADLIDARSR532_14)) then
            CFADLIDARSR532_14 = reshape(lidar_cfad_sr(:, 14, LM:1:-1), (/IM, JM, LM /))
         end if

         if (associated(CFADLIDARSR532_15)) then
            CFADLIDARSR532_15 = reshape(lidar_cfad_sr(:, 15, LM:1:-1), (/IM, JM, LM /))
         end if

         if (associated(LIDARPMOL)) then
            LIDARPMOL = reshape(lidar_pmol(:, LM:1:-1), (/IM, JM, LM /))
         end if

         if (associated(LIDARPTOT)) then
            LIDARPTOT = reshape(sum(lidar_beta_tot(:, :, LM:1:-1), 2) / (1. * ncolumns), (/IM, JM, LM /))
         end if

         if (associated(LIDARTAUTOT)) then
            LIDARTAUTOT = reshape(sum(lidar_tau_tot(:, :, LM:1:-1), 2) / (1. * ncolumns), (/IM, JM, LM /))
         end if

         ! CLOUDSAT/radar exports
         !-------------------------------------------------
         !-------------------------------------------------

#define SATSIMTEMP
#ifdef SATSIMTEMP
         if (associated(RADARZETOT)) then
            radar_ze_tot_mask = (radar_ze_tot /= MAPL_UNDEF) ! mask out undefined columns
            ncvalid = count(radar_ze_tot_mask, dim=2) ! number of defined columns
            radar_ze_tot_max = maxval(radar_ze_tot, dim=2, mask=radar_ze_tot_mask) ! maximum among columns
            radar_ze_tot_tmp = spread(radar_ze_tot_max, dim=2, ncopies=ncolumns) ! max copied to all columns
            where (radar_ze_tot_mask) &
                 radar_ze_tot_tmp = 10.0**((radar_ze_tot - radar_ze_tot_tmp) / 10.0) ! convert to rel reflectivity
            where (ncvalid > 0)
               radar_ze_tot_mean = &
                    sum(radar_ze_tot_tmp, dim=2, mask=radar_ze_tot_mask) / real(ncvalid) ! mean rel reflectivity
               radar_ze_tot_mean = radar_ze_tot_max + 10.0 * log10(radar_ze_tot_mean) ! convert back to dBZ
               elsewhere
               radar_ze_tot_mean = MAPL_UNDEF
            end where
            RADARZETOT = reshape(radar_ze_tot_mean(:, LM:1:-1), (/IM, JM, LM/))
         end if
         !if (associated(RADARZETOT)) then
         ! radar_ze_tot_mask = (radar_ze_tot .ne. MAPL_UNDEF)
         ! ncvalid = count( radar_ze_tot_mask, 2 )
         ! where (ncvalid > 0)
         !    radar_ze_tot_max = maxval( radar_ze_tot, 2, radar_ze_tot_mask )
         !    radar_ze_tot_mean = radar_ze_tot_max + 10.0*log10( sum( 10.0**((radar_ze_tot - spread(radar_ze_tot_max, 2,
         ! NCOLUMNS))/10.0), 2, radar_ze_tot_mask ) / ncvalid )
         ! elsewhere
         !    radar_ze_tot_mean = MAPL_UNDEF
         ! end where
         ! RADARZETOT=reshape( radar_ze_tot_mean(:,LM:1:-1), (/IM , JM , LM /) )
         !end if
#else
         if (associated(RADARZETOT)) then
            RADARZETOT = reshape(sum(radar_ze_tot(:, :, LM:1:-1), 2) / (1. * ncolumns), (/IM, JM, LM /))
         end if
#endif ! SATSIMTEMP

         if (associated(RADARLTCC)) then
            RADARLTCC = reshape(radar_lidar_tcc, (/IM, JM /))
         end if

         if (associated(CLCALIPSO2)) then
            CLCALIPSO2 = reshape(radar_lidar_only_freq_cloud(:, LM:1:-1), (/IM, JM, LM /))
         end if

         if (associated(CLOUDSATCFAD01)) then
            CLOUDSATCFAD01 = reshape(radar_cfad_ze(:, 1, LM:1:-1), (/IM, JM, LM /))
         end if

         if (associated(CLOUDSATCFAD02)) then
            CLOUDSATCFAD02 = reshape(radar_cfad_ze(:, 2, LM:1:-1), (/IM, JM, LM /))
         end if

         if (associated(CLOUDSATCFAD03)) then
            CLOUDSATCFAD03 = reshape(radar_cfad_ze(:, 3, LM:1:-1), (/IM, JM, LM /))
         end if

         if (associated(CLOUDSATCFAD04)) then
            CLOUDSATCFAD04 = reshape(radar_cfad_ze(:, 4, LM:1:-1), (/IM, JM, LM /))
         end if

         if (associated(CLOUDSATCFAD05)) then
            CLOUDSATCFAD05 = reshape(radar_cfad_ze(:, 5, LM:1:-1), (/IM, JM, LM /))
         end if

         if (associated(CLOUDSATCFAD06)) then
            CLOUDSATCFAD06 = reshape(radar_cfad_ze(:, 6, LM:1:-1), (/IM, JM, LM /))
         end if

         if (associated(CLOUDSATCFAD07)) then
            CLOUDSATCFAD07 = reshape(radar_cfad_ze(:, 7, LM:1:-1), (/IM, JM, LM /))
         end if

         if (associated(CLOUDSATCFAD08)) then
            CLOUDSATCFAD08 = reshape(radar_cfad_ze(:, 8, LM:1:-1), (/IM, JM, LM /))
         end if

         if (associated(CLOUDSATCFAD09)) then
            CLOUDSATCFAD09 = reshape(radar_cfad_ze(:, 9, LM:1:-1), (/IM, JM, LM /))
         end if

         if (associated(CLOUDSATCFAD10)) then
            CLOUDSATCFAD10 = reshape(radar_cfad_ze(:, 10, LM:1:-1), (/IM, JM, LM /))
         end if

         if (associated(CLOUDSATCFAD11)) then
            CLOUDSATCFAD11 = reshape(radar_cfad_ze(:, 11, LM:1:-1), (/IM, JM, LM /))
         end if

         if (associated(CLOUDSATCFAD12)) then
            CLOUDSATCFAD12 = reshape(radar_cfad_ze(:, 12, LM:1:-1), (/IM, JM, LM /))
         end if

         if (associated(CLOUDSATCFAD13)) then
            CLOUDSATCFAD13 = reshape(radar_cfad_ze(:, 13, LM:1:-1), (/IM, JM, LM /))
         end if

         if (associated(CLOUDSATCFAD14)) then
            CLOUDSATCFAD14 = reshape(radar_cfad_ze(:, 14, LM:1:-1), (/IM, JM, LM /))
         end if

         if (associated(CLOUDSATCFAD15)) then
            CLOUDSATCFAD15 = reshape(radar_cfad_ze(:, 15, LM:1:-1), (/IM, JM, LM /))
         end if

         ! MODIS exports -- currently available only from COSP

         if (associated(MDSCLDFRCTTL)) then
            MDSCLDFRCTTL = reshape(MODIS_Cloud_Fraction_Total_Mean, (/IM, JM /))
         end if

         if (associated(MDSCLDFRCWTR)) then
            MDSCLDFRCWTR = reshape(MODIS_Cloud_Fraction_Water_Mean, (/IM, JM /))
         end if

         if (associated(MDSCLDFRCH2O)) then
            MDSCLDFRCH2O = reshape(MODIS_Cloud_Fraction_Water_Mean, (/IM, JM /))
         end if

         if (associated(MDSCLDFRCICE)) then
            MDSCLDFRCICE = reshape(MODIS_Cloud_Fraction_Ice_Mean, (/IM, JM /))
         end if

         if (associated(MDSCLDFRCHI)) then
            MDSCLDFRCHI = reshape(MODIS_Cloud_Fraction_High_Mean, (/IM, JM /))
         end if

         if (associated(MDSCLDFRCMID)) then
            MDSCLDFRCMID = reshape(MODIS_Cloud_Fraction_Mid_Mean, (/IM, JM /))
         end if

         if (associated(MDSCLDFRCLO)) then
            MDSCLDFRCLO = reshape(MODIS_Cloud_Fraction_Low_Mean, (/IM, JM /))
         end if

         if (associated(MDSOPTHCKTTL)) then
            MDSOPTHCKTTL = reshape(MODIS_Optical_Thickness_Total_Mean, (/IM, JM /))
         end if

         if (associated(MDSOPTHCKWTR)) then
            MDSOPTHCKWTR = reshape(MODIS_Optical_Thickness_Water_Mean, (/IM, JM /))
         end if

         if (associated(MDSOPTHCKH2O)) then
            MDSOPTHCKH2O = reshape(MODIS_Optical_Thickness_Water_Mean, (/IM, JM /))
         end if

         if (associated(MDSOPTHCKICE)) then
            MDSOPTHCKICE = reshape(MODIS_Optical_Thickness_Ice_Mean, (/IM, JM /))
         end if

         if (associated(MDSOPTHCKTTLLG)) then
            MDSOPTHCKTTLLG = reshape(MODIS_Optical_Thickness_Total_LogMean, (/IM, JM /))
         end if

         if (associated(MDSOPTHCKWTRLG)) then
            MDSOPTHCKWTRLG = reshape(MODIS_Optical_Thickness_Water_LogMean, (/IM, JM /))
         end if

         if (associated(MDSOPTHCKH2OLG)) then
            MDSOPTHCKH2OLG = reshape(MODIS_Optical_Thickness_Water_LogMean, (/IM, JM /))
         end if

         if (associated(MDSOPTHCKICELG)) then
            MDSOPTHCKICELG = reshape(MODIS_Optical_Thickness_Ice_LogMean, (/IM, JM /))
         end if

         if (associated(MDSCLDSZWTR)) then
            MDSCLDSZWTR = reshape(MODIS_Cloud_Particle_Size_Water_Mean, (/IM, JM /))
         end if

         if (associated(MDSCLDSZH20)) then
            MDSCLDSZH20 = reshape(MODIS_Cloud_Particle_Size_Water_Mean, (/IM, JM /))
         end if

         if (associated(MDSCLDSZICE)) then
            MDSCLDSZICE = reshape(MODIS_Cloud_Particle_Size_Ice_Mean, (/IM, JM /))
         end if

         if (associated(MDSCLDTOPPS)) then
            MDSCLDTOPPS = reshape(MODIS_Cloud_Top_Pressure_Total_Mean, (/IM, JM /))
         end if

         if (associated(MDSWTRPATH)) then
            MDSWTRPATH = reshape(MODIS_Liquid_Water_Path_Mean, (/IM, JM /))
         end if

         if (associated(MDSH2OPATH)) then
            MDSH2OPATH = reshape(MODIS_Liquid_Water_Path_Mean, (/IM, JM /))
         end if

         if (associated(MDSICEPATH)) then
            MDSICEPATH = reshape(MODIS_Ice_Water_Path_Mean, (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST11)) then
            MDSTAUPRSHIST11 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 1, 1), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST12)) then
            MDSTAUPRSHIST12 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 1, 2), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST13)) then
            MDSTAUPRSHIST13 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 1, 3), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST14)) then
            MDSTAUPRSHIST14 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 1, 4), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST15)) then
            MDSTAUPRSHIST15 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 1, 5), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST16)) then
            MDSTAUPRSHIST16 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 1, 6), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST17)) then
            MDSTAUPRSHIST17 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 1, 7), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST21)) then
            MDSTAUPRSHIST21 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 2, 1), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST22)) then
            MDSTAUPRSHIST22 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 2, 2), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST23)) then
            MDSTAUPRSHIST23 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 2, 3), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST24)) then
            MDSTAUPRSHIST24 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 2, 4), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST25)) then
            MDSTAUPRSHIST25 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 2, 5), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST26)) then
            MDSTAUPRSHIST26 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 2, 6), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST27)) then
            MDSTAUPRSHIST27 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 2, 7), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST31)) then
            MDSTAUPRSHIST31 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 3, 1), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST32)) then
            MDSTAUPRSHIST32 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 3, 2), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST33)) then
            MDSTAUPRSHIST33 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 3, 3), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST34)) then
            MDSTAUPRSHIST34 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 3, 4), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST35)) then
            MDSTAUPRSHIST35 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 3, 5), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST36)) then
            MDSTAUPRSHIST36 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 3, 6), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST37)) then
            MDSTAUPRSHIST37 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 3, 7), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST41)) then
            MDSTAUPRSHIST41 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 4, 1), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST42)) then
            MDSTAUPRSHIST42 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 4, 2), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST43)) then
            MDSTAUPRSHIST43 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 4, 3), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST44)) then
            MDSTAUPRSHIST44 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 4, 4), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST45)) then
            MDSTAUPRSHIST45 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 4, 5), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST46)) then
            MDSTAUPRSHIST46 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 4, 6), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST47)) then
            MDSTAUPRSHIST47 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 4, 7), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST51)) then
            MDSTAUPRSHIST51 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 5, 1), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST52)) then
            MDSTAUPRSHIST52 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 5, 2), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST53)) then
            MDSTAUPRSHIST53 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 5, 3), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST54)) then
            MDSTAUPRSHIST54 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 5, 4), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST55)) then
            MDSTAUPRSHIST55 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 5, 5), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST56)) then
            MDSTAUPRSHIST56 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 5, 6), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST57)) then
            MDSTAUPRSHIST57 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 5, 7), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST61)) then
            MDSTAUPRSHIST61 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 6, 1), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST62)) then
            MDSTAUPRSHIST62 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 6, 2), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST63)) then
            MDSTAUPRSHIST63 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 6, 3), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST64)) then
            MDSTAUPRSHIST64 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 6, 4), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST65)) then
            MDSTAUPRSHIST65 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 6, 5), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST66)) then
            MDSTAUPRSHIST66 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 6, 6), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST67)) then
            MDSTAUPRSHIST67 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 6, 7), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST71)) then
            MDSTAUPRSHIST71 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 7, 1), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST72)) then
            MDSTAUPRSHIST72 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 7, 2), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST73)) then
            MDSTAUPRSHIST73 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 7, 3), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST74)) then
            MDSTAUPRSHIST74 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 7, 4), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST75)) then
            MDSTAUPRSHIST75 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 7, 5), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST76)) then
            MDSTAUPRSHIST76 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 7, 6), (/IM, JM /))
         end if

         if (associated(MDSTAUPRSHIST77)) then
            MDSTAUPRSHIST77 = reshape(MODIS_Optical_Thickness_vs_Cloud_Top_Pressure(:, 7, 7), (/IM, JM /))
         end if

         if (associated(MISRMNCLDTP)) then
            MISRMNCLDTP = reshape(MISR_meanztop, (/IM, JM /))
         end if

         if (associated(MISRCLDAREA)) then
            MISRCLDAREA = reshape(MISR_cldarea, (/IM, JM /))
         end if

         if (associated(MISRLYRTP0)) then
            MISRLYRTP0 = reshape(MISR_dist_model_layertops(:, 1), (/IM, JM /))
         end if

         if (associated(MISRLYRTP250)) then
            MISRLYRTP250 = reshape(MISR_dist_model_layertops(:, 2), (/IM, JM /))
         end if

         if (associated(MISRLYRTP750)) then
            MISRLYRTP750 = reshape(MISR_dist_model_layertops(:, 3), (/IM, JM /))
         end if

         if (associated(MISRLYRTP1250)) then
            MISRLYRTP1250 = reshape(MISR_dist_model_layertops(:, 4), (/IM, JM /))
         end if

         if (associated(MISRLYRTP1750)) then
            MISRLYRTP1750 = reshape(MISR_dist_model_layertops(:, 5), (/IM, JM /))
         end if

         if (associated(MISRLYRTP2250)) then
            MISRLYRTP2250 = reshape(MISR_dist_model_layertops(:, 6), (/IM, JM /))
         end if

         if (associated(MISRLYRTP2750)) then
            MISRLYRTP2750 = reshape(MISR_dist_model_layertops(:, 7), (/IM, JM /))
         end if

         if (associated(MISRLYRTP3500)) then
            MISRLYRTP3500 = reshape(MISR_dist_model_layertops(:, 8), (/IM, JM /))
         end if

         if (associated(MISRLYRTP4500)) then
            MISRLYRTP4500 = reshape(MISR_dist_model_layertops(:, 9), (/IM, JM /))
         end if

         if (associated(MISRLYRTP6000)) then
            MISRLYRTP6000 = reshape(MISR_dist_model_layertops(:, 10), (/IM, JM /))
         end if

         if (associated(MISRLYRTP8000)) then
            MISRLYRTP8000 = reshape(MISR_dist_model_layertops(:, 11), (/IM, JM /))
         end if

         if (associated(MISRLYRTP10000)) then
            MISRLYRTP10000 = reshape(MISR_dist_model_layertops(:, 12), (/IM, JM /))
         end if

         if (associated(MISRLYRTP12000)) then
            MISRLYRTP12000 = reshape(MISR_dist_model_layertops(:, 13), (/IM, JM /))
         end if

         if (associated(MISRLYRTP14000)) then
            MISRLYRTP14000 = reshape(MISR_dist_model_layertops(:, 14), (/IM, JM /))
         end if

         if (associated(MISRLYRTP16000)) then
            MISRLYRTP16000 = reshape(MISR_dist_model_layertops(:, 15), (/IM, JM /))
         end if

         if (associated(MISRLYRTP18000)) then
            MISRLYRTP18000 = reshape(MISR_dist_model_layertops(:, 16), (/IM, JM /))
         end if

         if (associated(MISRFQ0)) then
            MISRFQ0 = reshape(fq_MISR(:, :, 1), (/IM, JM, 7 /))
         end if

         if (associated(MISRFQ250)) then
            MISRFQ250 = reshape(fq_MISR(:, :, 2), (/IM, JM, 7 /))
         end if

         if (associated(MISRFQ750)) then
            MISRFQ750 = reshape(fq_MISR(:, :, 3), (/IM, JM, 7 /))
         end if

         if (associated(MISRFQ1250)) then
            MISRFQ1250 = reshape(fq_MISR(:, :, 4), (/IM, JM, 7 /))
         end if

         if (associated(MISRFQ1750)) then
            MISRFQ1750 = reshape(fq_MISR(:, :, 5), (/IM, JM, 7 /))
         end if

         if (associated(MISRFQ2250)) then
            MISRFQ2250 = reshape(fq_MISR(:, :, 6), (/IM, JM, 7 /))
         end if

         if (associated(MISRFQ2750)) then
            MISRFQ2750 = reshape(fq_MISR(:, :, 7), (/IM, JM, 7 /))
         end if

         if (associated(MISRFQ3500)) then
            MISRFQ3500 = reshape(fq_MISR(:, :, 8), (/IM, JM, 7 /))
         end if

         if (associated(MISRFQ4500)) then
            MISRFQ4500 = reshape(fq_MISR(:, :, 9), (/IM, JM, 7 /))
         end if

         if (associated(MISRFQ6000)) then
            MISRFQ6000 = reshape(fq_MISR(:, :, 10), (/IM, JM, 7 /))
         end if

         if (associated(MISRFQ8000)) then
            MISRFQ8000 = reshape(fq_MISR(:, :, 11), (/IM, JM, 7 /))
         end if

         if (associated(MISRFQ10000)) then
            MISRFQ10000 = reshape(fq_MISR(:, :, 12), (/IM, JM, 7 /))
         end if

         if (associated(MISRFQ12000)) then
            MISRFQ12000 = reshape(fq_MISR(:, :, 13), (/IM, JM, 7 /))
         end if

         if (associated(MISRFQ14000)) then
            MISRFQ14000 = reshape(fq_MISR(:, :, 14), (/IM, JM, 7 /))
         end if

         if (associated(MISRFQ16000)) then
            MISRFQ16000 = reshape(fq_MISR(:, :, 15), (/IM, JM, 7 /))
         end if

         if (associated(MISRFQ18000)) then
            MISRFQ18000 = reshape(fq_MISR(:, :, 16), (/IM, JM, 7 /))
         end if

         call ESMF_UserCompGetInternalState(GC, 'SatSim_State', wrap, status)
         _VERIFY(status)
         self => wrap%PTR

         if (self%nmask_vars > 0) then
            ! get orb bundle
            call ESMF_StateGet(IMPORT, 'SATORB', bundle, _RC)
            call MAPL_Get(STATE, ExportSpec=ExportSpec, _RC)
            do i = 1, self%nmask_vars
               vindex = MAPL_VarSpecGetIndex(ExportSpec, self%export_name(i), _RC)
               call MAPL_VarSpecGet(ExportSpec(vindex), DIMS=mapl_dims, ungridded_dims=ungridded_dims, _RC)
               ! mask from bundle
               call ESMFL_BundleGetPointerToData(bundle, self%mask_name(i), ptr_mask, _RC)
               if (self%newvar(i)) then

                  ! determine what dimension pointer we need
                  if (mapl_dims == MAPL_DimsHorzOnly .and. (.not. associated(ungridded_dims))) then
                     call MAPL_GetPointer(EXPORT, ptr2d_new, self%newvar_name(i), _RC)
                     if (associated(ptr2d_new)) then
                        call MAPL_GetPointer(EXPORT, ptr2d, self%export_name(i), _RC)
                        ptr2d_new = ptr2d
                        where (ptr_mask == MAPL_UNDEF)
                           ptr2d_new = MAPL_UNDEF
                        end where
                        nullify(ptr2d_new)
                        nullify(ptr2d)
                     end if
                  end if
                  if (mapl_dims == MAPL_DimsHorzOnly .and. associated(ungridded_dims)) then
                     call MAPL_GetPointer(EXPORT, ptr3d_new, self%newvar_name(i), _RC)
                     if (associated(ptr3d_new)) then
                        call MAPL_GetPointer(EXPORT, ptr3d, self%export_name(i), _RC)
                        ptr3d_new = ptr3d
                        do j = 1, ungridded_dims(1)
                           where (ptr_mask == MAPL_UNDEF)
                              ptr3d_new(:, :, j) = MAPL_UNDEF
                           end where
                           nullify(ptr3d_new)
                           nullify(ptr3d)
                        end do
                     end if
                  end if
                  if (mapl_dims == MAPL_DimsHorzVert) then
                     call MAPL_GetPointer(EXPORT, ptr3d_new, self%newvar_name(i), _RC)
                     if (associated(ptr3d_new)) then
                        call MAPL_GetPointer(EXPORT, ptr3d, self%export_name(i), _RC)
                        ptr3d_new = ptr3d
                        do j = 1, LM
                           where (ptr_mask == MAPL_UNDEF)
                              ptr3d_new(:, :, j) = MAPL_UNDEF
                           end where
                        end do
                        nullify(ptr3d_new)
                        nullify(ptr3d)
                     end if
                  end if

               else
                  ! determine what dimension pointer we need
                  if (mapl_dims == MAPL_DimsHorzOnly .and. (.not. associated(ungridded_dims))) then
                     call MAPL_GetPointer(EXPORT, ptr2d, self%newvar_name(i), _RC)
                     if (associated(ptr2d)) then
                        where (ptr_mask == MAPL_UNDEF)
                           ptr2d = MAPL_UNDEF
                        end where
                        nullify(ptr2d)
                     end if
                  end if
                  if (mapl_dims == MAPL_DimsHorzOnly .and. associated(ungridded_dims)) then
                     call MAPL_GetPointer(EXPORT, ptr3d, self%newvar_name(i), _RC)
                     if (associated(ptr3d)) then
                        do j = 1, ungridded_dims(1)
                           where (ptr_mask == MAPL_UNDEF)
                              ptr3d(:, :, j) = MAPL_UNDEF
                           end where
                        end do
                        nullify(ptr3d)
                     end if
                  end if
                  if (mapl_dims == MAPL_DimsHorzVert) then
                     call MAPL_GetPointer(EXPORT, ptr3d, self%newvar_name(i), _RC)
                     if (associated(ptr3d)) then
                        do j = 1, LM
                           where (ptr_mask == MAPL_UNDEF)
                              ptr3d(:, :, j) = MAPL_UNDEF
                           end where
                        end do
                        nullify(ptr3d)
                     end if
                  end if

               end if

               nullify(ptr_mask)

            end do
         end if

         ! deallocate stuff
         !-----------------
         deallocate(frac_out, _STAT)
         !deallocate(      frac_outinv, _STAT)
         deallocate(lidar_beta_tot, _STAT)
         deallocate(lidar_tau_tot, _STAT)
         deallocate(radar_ze_tot, _STAT)
         deallocate(radar_ze_tot_tmp, _STAT)
         deallocate(radar_ze_tot_mask, _STAT)

         call FREE_COSP_GRIDBOX(gbx)
         call FREE_COSP_SUBGRID(sgx)
         if (cfg%Lisccp_sim) call FREE_COSP_ISCCP(isccp)
         !call FREE_COSP_SGHYDRO(sghydro)
         if (cfg%Lradar_sim) call FREE_COSP_SGRADAR(sgradar)
         if (cfg%Llidar_sim) call FREE_COSP_SGLIDAR(sglidar)
         call FREE_COSP_VGRID(vgrid)
         if (cfg%Llidar_sim) call FREE_COSP_LIDARSTATS(stlidar)
         if (cfg%Lradar_sim) call FREE_COSP_RADARSTATS(stradar)
         if (cfg%Lmodis_sim) call FREE_COSP_MODIS(modis)
         if (cfg%Lmisr_sim) call FREE_COSP_MISR(misr)

         _RETURN(_SUCCESS)

      end subroutine SIM_DRIVER

   end subroutine Run

end module GEOS_SatsimGridCompMod
