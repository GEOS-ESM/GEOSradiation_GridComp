#include "MAPL_Generic.h"

module GEOS_IrradGridCompMod

   !BOP
   !MODULE: GEOS_Irrad -- Computes longwave radiative transfer through a cloudy atmosphere

   !DESCRIPTION:
   !
   !   {\tt Irrad} is a light-weight gridded component to compute longwave
   ! radiative fluxes. It operates on the ESMF grid that appears in the
   ! gridded component passed to its {\tt Initialize} method. Unlike
   ! heavier gridded components, it does not enforce its own grid.
   ! The only restrictions are that it be a 3-dimensional grid
   ! in which one dimension is aligned with the vertical coordinate and
   ! only the horizontal dimensions are decomposed.
   !
   !   The radiative transfer calculation is based on M-D Chou's IRRAD routine.
   ! A full documentation of the code may be found in
   ! "A Thermal Infrared Radiation Parameterization for Atmospheric Studies"
   ! M.-D. Chou et al., NASA/TM-2001-104606, Vol. 19, 55 pp, 2003.
   ! Based on the 1996-version of the Air Force Geophysical Laboratory HITRAN data
   ! base (Rothman et al., 1998), the parameterization includes the absorption due
   ! to major gaseous absorption (water vapor, CO2 , O3 ) and most of the minor
   ! trace gases (N2O, CH4 , CFC's), as well as clouds and aerosols. The thermal
   ! infrared spectrum is divided into nine bands and a subband. To achieve a high
   ! degree of accuracy and speed, various approaches of computing the transmission
   ! function are applied to different spectral bands and gases. The gaseous
   ! transmission function is computed either using the k-distribution method or
   ! the table look-up method. To include the effect of scattering due to clouds
   ! and aerosols, the optical thickness is scaled by the single-scattering albedo
   ! and asymmetry factor. The optical thickness, the single-scattering albedo,
   ! and the asymmetry factor of clouds are parameterized as functions of the ice
   ! and water content and the particle size.

   !   All outputs are optional and are filled only if they have been
   ! initialized by a coupler.
   !
   !   The net (+ve downward) fluxes are returned at the layer
   ! interfaces, which are indexed from the top of the atmosphere (L=0)
   ! to the surface. It also computes the sensitivity of net downward flux to
   ! surface temperature and emission by the surface.
   ! The full transfer calculation, including the linearization w.r.t. the surface temperature,
   ! is done intermitently, on the component's main time step and its results are
   ! kept in the internal state. Exports are refreshed each heartbeat based on the
   ! latest surface temperature.
   !
   !   Radiation should be called either before or after thos components
   !    (usually SURFACE and DYNAMICS) that use its fluxes and modify
   !    its inputs. If it is called before, the intemittent refresh should
   !    occur during the first step of the radiation cycle, while if it
   !    is called after, it should occur during the last step. The behavior
   !    of the component needs to be somewhat different in these two cases
   !    and so a means is provided, through the logical attribute \texttt{CALL\_LAST} in
   !    configuration, of telling the component how it is being used. The
   !    default is \texttt{CALL\_LAST = "TRUE"}.

   !USES:
   use ESMF
   use MAPL
   use GEOS_UtilsMod
   use gFTL_StringVector

   use rrtmg_lw_rad, only: rrtmg_lw
   use rrtmg_lw_init, only: rrtmg_lw_ini
   use parrrtm, only: ngptlw, nbndlw
   use rrlw_wvn, only: wavenum1, wavenum2
   use rad_utils, only: Tbr_from_band_flux, choose_solar_scheme, choose_irrad_scheme

   ! for RRTMGP
   use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp

   use irradmod, only: IRRAD

   implicit none
   private

   !PUBLIC MEMBER FUNCTIONS:

   public SetServices

   !EOP

   ! number of bands in different radiation codes

   integer, parameter :: NB_CHOU   = 10        ! #bands in IRRAD calcs for Chou
   integer, parameter :: NB_RRTMG  = 16        ! #bands in IRRAD calcs for RRTMG
   integer, parameter :: NB_RRTMGP = 16        ! #bands in IRRAD calcs for RRTMGP

   integer, parameter :: NB_CHOU_SORAD   = 8   ! #bands in SORAD calcs for Chou
   integer, parameter :: NB_RRTMG_SORAD  = 14  ! #bands in SORAD calcs for RRTMG
   integer, parameter :: NB_RRTMGP_SORAD = 14  ! #bands in SORAD calcs for RRTMGP

   ! Select which RRTMG[P] bands support OLR output
   !    via OLRBbbRG and TBRBbbRG exports ...
   ! (These exports require support space in the
   ! internal state so we choose only the ones we want
   ! to offer here at compile time. In the future, if
   ! most of the bands are required, then we will change
   ! strategy and reserve space for ALL of them in the
   ! internal state. In that case we would only require
   ! runtime band selection via the EXPORTS chosen.)

   ! NOTE: RRTMG band 16 should be requested with caution ...
   ! Band 16 is technically 2600-3250 cm-1. But when RT
   ! is performed across all 16 bands, as it is in GEOS-5
   ! usage, then band 16 includes the integrated Planck
   ! values from 2600 cm-1 to infinity. So, the brightness
   ! temperature (Tbr) calculations in Update_Flx(), which
   ! use specified wavenumber endpoints, may require mod-
   ! ification for band 16 (pmn: TODO). For the moment,
   ! the limits [2600,3250] are used.

   ! Which bands are supported?
   !    (Currently RRTMG & RRTMGP only:
   !      RRTMG & RRTMGP have the same number of bands &
   !      very similar, but not identical, band limits)
   !    (Actual calculation only if export is requested)
   ! Supported?    Band  Requested by (and use)
   logical, parameter :: band_output_supported (nbndlw) = [ &
        .false. , &!  01
        .false. , &!  02
        .false. , &!  03
        .false. , &!  04
        .true.  , &!  05   W. Putman (CO2 Longwave IR, GOES Band 16)
        .true.  , &!  06   A. Collow (Longwave IR, GOES Band 14)
        .true.  , &!  07   W. Putman (Ozone IR, GOES Band 12)
        .true.  , &!  08   W. Putman (needed for lightning param)
        .true.  , &!  09   W. Putman (Lower-level Water Vapor, GOES Band 10)
        .true.  , &!  10   W. Putman (Mid-level Water Vapor, GOES Band 9)
        .true.  , &!  11   W. Putman (Upper-level Water Vapor, GOES Band 8)
        .false. , &!  12
        .false. , &!  13
        .false. , &!  14
        .true.  , &!  15   W. Putman (Shortwave IR, GOES Band 7)
        .false. ]  !  16
   ! PMN: TODO, make LW method like SW so it doesnt waste
   ! intermediate variable space on unused bands?

   ! PS: We may later have an RRTMG internal state like
   ! RRTMGP below with various rrtmg_lw_init data, etc.
   ! TODO

   ! RRTMGP internal state
   ! This will be attached to the Gridded Component
   ! used to provide efficient initialization
   type ty_RRTMGP_state
      private
      logical :: initialized = .false.
      type (ty_gas_optics_rrtmgp) :: k_dist
   end type ty_RRTMGP_state

   ! wrapper to access RRTMGP internal state
   type ty_RRTMGP_wrap
      type (ty_RRTMGP_state), pointer :: ptr => null()
   end type ty_RRTMGP_wrap

contains

   !BOP
   !IROUTINE: SetServices -- Sets ESMF services for this component

   !INTERFACE:
   subroutine SetServices ( GC, RC )

      !ARGUMENTS:
      type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
      integer, optional                  :: RC  ! return code

      !DESCRIPTION: This version uses the MAPL\_GenericSetServices. This function sets
      !                the Initialize and Finalize services, as well as allocating
      !   our instance of a generic state and putting it in the
      !   gridded component (GC). Here we only need to set the run method and
      !   add the state variable specifications (also generic) to our instance
      !   of the generic state. This is the way our true state variables get into
      !   the ESMF\_State INTERNAL, which is in the MAPL\_MetaComp.

      !EOP

      character(len=ESMF_MAXSTR) :: IAm
      integer                    :: STATUS
      character(len=ESMF_MAXSTR) :: COMP_NAME

      type (MAPL_MetaComp), pointer :: MAPL
      type (ESMF_Config) :: CF

      integer :: MY_STEP
      integer :: ACCUMINT
      real    :: DT

      logical :: USE_RRTMGP, USE_RRTMG, USE_CHOU

      type (ty_RRTMGP_state), pointer :: rrtmgp_state
      type (ty_RRTMGP_wrap)           :: wrap

      ! <<>> MSL
      integer      :: i,n
      character(len=ESMF_MAXSTR) :: gen_str !i.e. generic_string variable <<>> MSL
      character(len=ESMF_MAXSTR), allocatable :: nameRATS(:)

      ! Get my name and set-up traceback handle
      call ESMF_GridCompGet(GC, NAME=COMP_NAME, _RC)
      Iam = trim(COMP_NAME) // 'SetServices'

      ! save pointer to the wrapped RRTMGP internal state in the GC
      allocate(rrtmgp_state, _STAT)
      wrap%ptr => rrtmgp_state
      call ESMF_UserCompSetInternalState(GC, 'RRTMGP_state', wrap, status)
      VERIFY_(status)

      ! Get my internal MAPL_Generic state
      call MAPL_GetObjectFromGC (GC, MAPL, _RC)

      ! Get the intervals; "heartbeat" must exist
      call MAPL_GetResource (MAPL, DT, Label="RUN_DT:", _RC)

      ! Refresh interval defaults to heartbeat.
      ! This will also be read by MAPL_Generic and set as the component's main time step.
      call MAPL_GetResource (MAPL, DT, Label=trim(COMP_NAME)//"_DT:", default=DT, _RC)
      MY_STEP = nint(DT)

      ! Averaging interval defaults to refresh interval.
      call MAPL_GetResource (MAPL, DT, Label=trim(COMP_NAME)//"Avrg:", default=DT, _RC)
      ACCUMINT = nint(DT)

      ! Decide which radiation to use
      call choose_irrad_scheme (MAPL, USE_RRTMGP, USE_RRTMG, USE_CHOU, _RC)

      ! Set the state variable specs.

      !BOS

      !IMPORT STATE:

#include "Irrad_Import___.h"

      call MAPL_AddImportSpec(GC,                                   &
           LONG_NAME  = 'aerosols',                                   &
           UNITS      = 'kg kg-1',                                    &
           SHORT_NAME = 'AERO',                                       &
           DIMS       = MAPL_DimsHorzVert,                            &
           VLOCATION  = MAPL_VLocationCenter,                         &
           DATATYPE   = MAPL_StateItem,                               &
           RESTART    = MAPL_RestartSkip,                      _RC )

      ! If CO2 is provided as a RAT, import a CO2 field <<>> MSL
      ! Using DT below since it is already declared, and avoids adding an additional var - MSL
      call ESMF_GridCompGet(GC, CONFIG=CF, _RC)
      call ESMF_ConfigGetAttribute(CF, DT, Label='CO2:', default=-1.0, _RC)

      ! If using 3-D CO2, set up the import
      if (DT.eq.-2.0) then
         call ESMF_ConfigFindLabel( CF,'CO2_PROVIDER',rc=RC )
         n = ESMF_ConfigGetLen(CF,label='CO2_PROVIDER',rc=status)
         call ESMF_ConfigFindLabel( CF,'CO2_PROVIDER',rc=RC ) ! Godda reset!
         call ESMF_ConfigGetAttribute(CF, gen_str, Label='CO2_PROVIDER', default='none', RC=STATUS)
         ! If CO2: is invalid, raise an error
         if (n .le. 0 .or. ESMF_UtilStringLowerCase(trim(gen_str)) .eq. 'none') then
            gen_str = 'In AGCM.rc, cannot set CO2: to -2 and not give a valid CO2_PROVIDER'
            __raise__(MAPL_RC_ERROR, gen_str)
         endif

         call MAPL_AddImportSpec(GC,                                  &
              SHORT_NAME         = 'CO2',                               &
              LONG_NAME          = 'carbondioxide_concentration',       &
              UNITS              = 'pppv',                              &
              DIMS               = MAPL_DimsHorzVert,                   &
              VLOCATION          = MAPL_VLocationCenter,                &
              AVERAGING_INTERVAL = ACCUMINT,                            &
              REFRESH_INTERVAL   = MY_STEP,                             &
              _RC)
      endif

      !EXPORT STATE:

#include "Irrad_Export___.h"

      !  Irrad does not have a "real" internal state. To update the net_longwave_flux
      !  due to the change of surface temperature every time step, we keep
      !  several variables in the internal state.

      !INTERNAL STATE:

#include "Irrad_Internal___.h"

      ! Settings for RATS-specific radiation diagnostics
      ! -- these will cause RRTMG_LW to be called multiple times toggling the named species
      !    on/off. Outputs the specific radiative impacts of that species (e.g. CO2)
      ! -- the code below simply adds the diagnostic exports
      ! -- Need to test if a given input has an import. Not sure how to do this within SetServices. In
      !    case it's possible, this code block is placed -after- all the imports are added so the
      !    imports can be queried.
      ! <<>> MSL
      call ESMF_ConfigFindLabel(CF, 'RATS_DIAGNOSTICS:', RC=STATUS) ! Use STATUS to test if label was found

      IF (STATUS .eq. ESMF_SUCCESS) THEN
         n = ESMF_ConfigGetLen(CF,label='RATS_DIAGNOSTICS:',_RC)
      ENDIF

      ! No error thrown. Just go around this if nothing learnable from config.
      IF (STATUS .eq. ESMF_SUCCESS .and. n .ne. 0 ) THEN ! if the label was found...

         ! Get number of words in config line
         n = ESMF_ConfigGetLen(CF,label='RATS_DIAGNOSTICS:',_RC)

         allocate(nameRATS(n), STAT=STATUS)
         VERIFY_(STATUS)

         ! Put the cursor at the label
         call ESMF_ConfigFindLabel(CF, 'RATS_DIAGNOSTICS:', _RC)

         ! Loop over RATS in list
         DO i=1,n
            call ESMF_ConfigGetAttribute(CF,gen_str,_RC)

            nameRATS(i) = trim(gen_str)

            ! Test if label 'gen_str' has an import
            ! Otherwise throw error. (This could be done in GEOS_Physics GC)
            ! <<TBD>> MSL
            ! O3 CO2 CH4 N2O CFC11 CFC12 CFC22 CCl4

         ENDDO

         DO i=1,n
            ! Can't read the CF list in the loop above if MAPL_AddExportSpec() is
            ! called within it. Not sure why. So...

            ! Create exports for this RAT
            ! -- OLR
            !          gen_str = 'OLR_'//trim(gen_str)
            call MAPL_AddExportSpec(GC,                                            &
                 SHORT_NAME = 'dOLR_'//trim(nameRATS(i)),                          &
                 LONG_NAME  = 'chg_in_upwell_LW_flx_at_toa_from_'//trim(nameRATS(i)), &
                 UNITS      = 'W m-2',                                             &
                 DIMS       = MAPL_DimsHorzOnly,                                   &
                 VLOCATION  = MAPL_VLocationNone,                       _RC)
            call MAPL_AddExportSpec(GC,                                    &
                 SHORT_NAME = 'dLWS_'//trim(nameRATS(i)),                  &
                 LONG_NAME  = 'chg_in_surface_absorbed_LW_rad_from_'//trim(nameRATS(i)), &
                 UNITS      = 'W m-2',                                     &
                 DIMS       = MAPL_DimsHorzOnly,                           &
                 VLOCATION  = MAPL_VLocationNone,               _RC)
            call MAPL_AddExportSpec(GC,                                    &
                 SHORT_NAME = 'dFLNS_'//trim(nameRATS(i)),                 &
                 LONG_NAME  = 'chg_in_sfc_net_downward_LW_flux_from_'//trim(nameRATS(i)),&
                 UNITS      = 'W m-2',                                     &
                 DIMS       = MAPL_DimsHorzOnly,                           &
                 VLOCATION  = MAPL_VLocationNone,               _RC)
            call MAPL_AddExportSpec(GC,                                  &
                 SHORT_NAME = 'dSFCEM_'//trim(nameRATS(i)),                &
                 LONG_NAME  = 'LW_flux_emitted_from_sfc_from_'//trim(nameRATS(i)), &
                 UNITS      = 'W m-2',                                     &
                 DIMS       = MAPL_DimsHorzOnly,                           &
                 VLOCATION  = MAPL_VLocationNone,               _RC)
            call MAPL_AddExportSpec(GC,                                    &
                 SHORT_NAME = 'NETTRAP_'//trim(nameRATS(i)),                 &
                 LONG_NAME  = 'Net_Heat_trapping_due_to_'//trim(nameRATS(i)),&
                 UNITS      = 'W m-2',                                     &
                 DIMS       = MAPL_DimsHorzOnly,                           &
                 VLOCATION  = MAPL_VLocationNone,               _RC)
            call MAPL_AddExportSpec(GC,                                    &
                 SHORT_NAME = 'COLTRAP_'//trim(nameRATS(i)),               &
                 LONG_NAME  = 'Heat_trapping_due_to_'//trim(nameRATS(i)),  &
                 UNITS      = 'W m-2',                                     &
                 DIMS       = MAPL_DimsHorzVert,                           &
                 VLOCATION  = MAPL_VLocationCenter,             _RC)

            call MAPL_AddExportSpec(GC,                                    &
                 SHORT_NAME = 'FLX_'//trim(nameRATS(i)),                   &
                 LONG_NAME  = 'net_downward_longwave_flux_in_air_due_to'//trim(nameRATS(i)), &
                 UNITS      = 'W m-2',                                     &
                 DIMS       = MAPL_DimsHorzVert,                           &
                 VLOCATION  = MAPL_VLocationEdge,                   _RC )

            call MAPL_AddInternalSpec(GC,                                  &
                 SHORT_NAME = 'DFDTS_'//trim(nameRATS(i)),                 &
                 LONG_NAME  = 'sensitivity_of_net_downward_longwave_flux_in_air_to_surface_temperature_due_to'//trim(nameRATS(i)),&
                 UNITS      = 'W m-2 K-1',                                 &
                 DIMS       = MAPL_DimsHorzVert,                           &
                 add2export = .true., &
                 VLOCATION  = MAPL_VLocationEdge,                   _RC )

         ENDDO
         call MAPL_AddExportSpec(GC,                                       &
              SHORT_NAME = 'CO2_FIXED',                                    &
              LONG_NAME  = 'lol',                                          &
              UNITS      = 'mol/mol',                                      &
              DIMS       = MAPL_DimsHorzVert,                              &
              VLOCATION  = MAPL_VLocationCenter,               _RC)
         call MAPL_AddExportSpec(GC,                                       &
              SHORT_NAME = 'DELT',                                         &
              LONG_NAME  = 'change in surface temperature in RRTMG',       &
              UNITS      = 'K',                                            &
              DIMS       = MAPL_DimsHorzOnly,                              &
              VLOCATION  = MAPL_VLocationCenter,               _RC)
         if (allocated(nameRATS)) deallocate(nameRATS, STAT=STATUS)
         VERIFY_(STATUS)
         ! end rats code <<>> MSL

         ! Add necessary internal fields
         call MAPL_AddInternalSpec(GC,                             &
              SHORT_NAME = 'FLXU_RAT',                                  &
              LONG_NAME  = 'upward_longwave_flux_in_air',               &
              UNITS      = 'W m-2',                                     &
              UNGRIDDED_DIMS     = (/N/),                               &
              DIMS       = MAPL_DimsHorzVert,                           &
              VLOCATION  = MAPL_VLocationEdge,                   _RC )

         call MAPL_AddInternalSpec(GC,                             &
              SHORT_NAME = 'FLXD_RAT',                                  &
              LONG_NAME  = 'upward_longwave_flux_in_air',               &
              UNITS      = 'W m-2',                                     &
              UNGRIDDED_DIMS     = (/N/),                               &
              DIMS       = MAPL_DimsHorzVert,                           &
              VLOCATION  = MAPL_VLocationEdge,                   _RC )

         call MAPL_AddInternalSpec(GC,                             &
              SHORT_NAME = 'FLX_RAT',                                   &
              LONG_NAME  = 'net_downward_longwave_flux_in_air',         &
              UNITS      = 'W m-2',                                     &
              UNGRIDDED_DIMS     = (/N/),                               &
              DIMS       = MAPL_DimsHorzVert,                           &
              VLOCATION  = MAPL_VLocationEdge,                   _RC )

         call MAPL_AddInternalSpec(GC,                             &
              SHORT_NAME = 'DFDTS_RAT',                                 &
              LONG_NAME  = 'sensitivity_of_net_downward_longwave_flux_in_air_to_surface_temperature', &
              UNITS      = 'W m-2 K-1',                                 &
              UNGRIDDED_DIMS     = (/N/),                               &
              DIMS       = MAPL_DimsHorzVert,                           &
              VLOCATION  = MAPL_VLocationEdge,                   _RC )

         call MAPL_AddInternalSpec(GC,                             &
              SHORT_NAME = 'SFCEM_RAT',                                 &
              LONG_NAME  = 'longwave_flux_emitted_from_surface',        &
              UNITS      = 'W m-2',                                     &
              UNGRIDDED_DIMS     = (/N/),                               &
              DIMS       = MAPL_DimsHorzOnly,                           &
              VLOCATION  = MAPL_VLocationNone,                   _RC )
      ENDIF

      !EOS

      ! Set the Profiling timers
      call MAPL_TimerAdd(GC, name="-LW_DRIVER"               , _RC)
      call MAPL_TimerAdd(GC, name="--IRRAD"                  , _RC)
      call MAPL_TimerAdd(GC, name="---IRRAD_RUN"             , _RC)
      call MAPL_TimerAdd(GC, name="--RRTMG"                  , _RC)
      call MAPL_TimerAdd(GC, name="---RRTMG_RUN"             , _RC)
      call MAPL_TimerAdd(GC, name="---RRTMG_INIT"            , _RC)
      call MAPL_TimerAdd(GC, name="---RRTMG_FLIP"            , _RC)
      call MAPL_TimerAdd(GC, name="--RRTMGP"                 , _RC)
      call MAPL_TimerAdd(GC, name="---RRTMGP_IO_GAS"         , _RC)
      call MAPL_TimerAdd(GC, name="---RRTMGP_IO_CLOUDS"      , _RC)
      call MAPL_TimerAdd(GC, name="---RRTMGP_CLOUD_OPTICS"   , _RC)
      call MAPL_TimerAdd(GC, name="---RRTMGP_MCICA"          , _RC)
      call MAPL_TimerAdd(GC, name="---RRTMGP_GAS_OPTICS"     , _RC)
      call MAPL_TimerAdd(GC, name="---RRTMGP_RT"             , _RC)
      call MAPL_TimerAdd(GC, name="---RRTMGP_POST"           , _RC)
      call MAPL_TimerAdd(GC, name="--MISC"                   , _RC)
      call MAPL_TimerAdd(GC, name="---AEROSOLS"              , _RC)
      call MAPL_TimerAdd(GC, name="-UPDATE_FLX"              , _RC)

      ! Set Run method and use generic Initalize and Finalize methods
      call MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_RUN, Run, _RC)
      call MAPL_GenericSetServices(GC, _RC)

      RETURN_(ESMF_SUCCESS)
   end subroutine SetServices

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !BOP
   !IROUTINE: RUN -- Run method for the LW component

   !INTERFACE:
   subroutine RUN ( GC, IMPORT, EXPORT, CLOCK, RC )

      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component
      type(ESMF_State),    intent(inout) :: IMPORT ! Import state
      type(ESMF_State),    intent(inout) :: EXPORT ! Export state
      type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
      integer, optional,   intent(  out) :: RC     ! Error code:

      !DESCRIPTION: Periodically refreshes the fluxes and their derivatives
      !                w.r.t surface skin temperature. On every step it produces
      !                a linear estimate of the fluxes based on the instantaneous
      !                surface temperature.

      !EOP

      ! ErrLog Variables

      character(len=ESMF_MAXSTR)         :: IAm
      integer                            :: STATUS
      character(len=ESMF_MAXSTR)         :: COMP_NAME

      ! Local derived type aliases

      type (MAPL_MetaComp),     pointer  :: MAPL
      type (ESMF_Grid)                   :: ESMFGRID
      type (ESMF_State)                  :: INTERNAL
      type (ESMF_Alarm)                  :: ALARM

      integer                            :: IM, JM, LM
      integer                            :: CalledLast

      type (ty_RRTMGP_state), pointer    :: rrtmgp_state => null()
      type (ty_RRTMGP_wrap)              :: wrap

      ! Pointers to internal

      real, pointer, dimension(:,:  )   :: SFCEM_INT
      real, pointer, dimension(:,:  )   :: TS_INT
      real, pointer, dimension(:,:,:)   :: FLX_INT
      real, pointer, dimension(:,:,:)   :: FLXA_INT
      real, pointer, dimension(:,:,:)   :: FLC_INT
      real, pointer, dimension(:,:,:)   :: FLA_INT
      real, pointer, dimension(:,:,:)   :: FLXU_INT
      real, pointer, dimension(:,:,:)   :: FLXAU_INT
      real, pointer, dimension(:,:,:)   :: FLCU_INT
      real, pointer, dimension(:,:,:)   :: FLAU_INT
      real, pointer, dimension(:,:,:)   :: FLXD_INT
      real, pointer, dimension(:,:,:)   :: FLXAD_INT
      real, pointer, dimension(:,:,:)   :: FLCD_INT
      real, pointer, dimension(:,:,:)   :: FLAD_INT

      real, pointer, dimension(:,:,:)   :: DFDTS
      real, pointer, dimension(:,:,:)   :: DFDTSNA
      real, pointer, dimension(:,:,:)   :: DFDTSC
      real, pointer, dimension(:,:,:)   :: DFDTSCNA

      real, external :: getco2

      ! Concerning what radiation to use (global to LW_driver and Update_Flx)

      logical :: USE_RRTMGP, USE_RRTMGP_SORAD
      logical :: USE_RRTMG,  USE_RRTMG_SORAD
      logical :: USE_CHOU,   USE_CHOU_SORAD

      integer :: NB_IRRAD  ! Number of bands in IRRAD calcs

      ! local
      integer :: TOTAL_RAD_BANDS, NUM_BANDS

      ! Additional pointers for RRTMG

      real, pointer, dimension(:,:  )   :: LONS
      real, pointer, dimension(:,:  )   :: LATS

      ! which bands require OLR output?
      ! (only RRTMG[P]; OLRBbbRG, TBRBbbRG)
      real, pointer, dimension(:,:) :: ptr2d
      logical :: band_output (nbndlw)
      logical :: any_band_output
      integer :: ibnd
      character*2 :: bb

      ! For RATS <<>> MSL
      type (ESMF_Config)              :: CF
      logical, save                   :: first = .true. ! I don't wanna do this, but there's no Initialize() method. This prevents repeating unneeded ops
      integer, save                   :: nRATS ! Number of active RATs to toggle
      character(len=6), dimension(8)  :: RATNAMES = (/'O3    ','N2O   ','CFC11 ','CFC12 ','CH4   ','HCFC22','H2O   ','CO2   '/)
      character(len=128)              :: gen_str
      character(len=128), allocatable, save :: nameRATS(:) ! Current toggles read from AGCM.rc
      real, allocatable, dimension(:,:)     :: TMP_R
      real, allocatable, dimension(:,:,:)   :: UFLXRAT, DFLXRAT, DUFLX_DT_RAT
      real,     pointer, dimension(:,:,:)   :: SFCEM_INT_RAT
      real,     pointer, dimension(:,:,:,:) :: DFDTS_RAT, FLX_INT_RAT, FLXU_INT_RAT, FLXD_INT_RAT

      ! Begin...

      ! Get the target components name and set-up traceback handle.

      Iam = "Run"
      call ESMF_GridCompGet( GC, name=COMP_NAME, GRID=ESMFGRID, CONFIG=CF, _RC)
      Iam = trim(COMP_NAME) // Iam

      ! Get my internal MAPL_Generic state

      call MAPL_GetObjectFromGC ( GC, MAPL, _RC)

      call MAPL_TimerOn(MAPL,"TOTAL")

      ! Get parameters from generic state. The RUNALARM is used to control
      !  the calling of the full transfer calculation

      call MAPL_Get(MAPL, &
           IM=IM, JM=JM, LM=LM, &
           LONS=LONS, LATS=LATS, &
           RUNALARM=ALARM, &
           INTERNAL_ESMF_STATE=INTERNAL, &
           _RC)

      ! Decide which radiation to use:
      ! These USE_ flags are shared globally by contained LW_Driver() and Update_Flx()
      call choose_irrad_scheme (MAPL, &
           USE_RRTMGP,       USE_RRTMG,       USE_CHOU,       _RC)
      call choose_solar_scheme (MAPL, &
           USE_RRTMGP_SORAD, USE_RRTMG_SORAD, USE_CHOU_SORAD, _RC)

      ! Set number of IRRAD bands
      if (USE_RRTMGP) then
         NB_IRRAD = NB_RRTMGP
      else if (USE_RRTMG) then
         NB_IRRAD = NB_RRTMG
      else
         NB_IRRAD = NB_CHOU
      end if

      ! Test to see if AGCM.rc is set up correctly for the Radiation selected
      TOTAL_RAD_BANDS = NB_IRRAD
      if (USE_RRTMGP_SORAD) then
         TOTAL_RAD_BANDS = TOTAL_RAD_BANDS + NB_RRTMGP_SORAD
      else if (USE_RRTMG_SORAD) then
         TOTAL_RAD_BANDS = TOTAL_RAD_BANDS + NB_RRTMG_SORAD
      else
         TOTAL_RAD_BANDS = TOTAL_RAD_BANDS + NB_CHOU_SORAD
      end if

      call MAPL_GetResource (MAPL, NUM_BANDS ,'NUM_BANDS:', _RC)
      if (NUM_BANDS /= TOTAL_RAD_BANDS) then
         if (MAPL_am_I_Root()) then
            write (*,*) "NUM_BANDS is not set up correctly for the radiation combination selected:"
            write (*,*) "    IRRAD RRTMG: ", USE_RRTMGP      , USE_RRTMG      , USE_CHOU
            write (*,*) "    SOLAR RRTMG: ", USE_RRTMGP_SORAD, USE_RRTMG_SORAD, USE_CHOU_SORAD
            write (*,*) "Please check that your optics tables and NUM_BANDS are correct."
         end if
         _FAIL('Total number of radiation bands is inconsistent!')
      end if

      ! select which bands require OLRB output ...
      ! Only available for RRTMG[P]
      ! must be supported AND requested by export 'OLRBbbRG' OR 'TBRBbbRG'
      any_band_output = .false.
      if (USE_RRTMG .or. USE_RRTMGP) then
         do ibnd = 1,nbndlw
            band_output(ibnd) = .false.
            if (.not. band_output_supported(ibnd)) cycle
            write(bb,'(I0.2)') ibnd
            call MAPL_GetPointer(EXPORT, ptr2d, 'OLRB'//bb//'RG', _RC)
            if (associated(ptr2d)) then
               band_output(ibnd) = .true.
               cycle
            end if
            call MAPL_GetPointer(EXPORT, ptr2d, 'TBRB'//bb//'RG', _RC)
            if (associated(ptr2d)) then
               band_output(ibnd) = .true.
               cycle
            end if
         end do
         any_band_output = any(band_output)
      end if

      ! Pointers to Internals; these are needed by both Update and Refresh

      call MAPL_GetPointer(INTERNAL, SFCEM_INT, 'SFCEM', _RC)
      call MAPL_GetPointer(INTERNAL, FLX_INT,   'FLX',   _RC)
      call MAPL_GetPointer(INTERNAL, FLXA_INT,  'FLXA',  _RC)
      call MAPL_GetPointer(INTERNAL, FLC_INT,   'FLC',   _RC)
      call MAPL_GetPointer(INTERNAL, FLA_INT,   'FLA',   _RC)
      call MAPL_GetPointer(INTERNAL, FLXU_INT,  'FLXU',  _RC)
      call MAPL_GetPointer(INTERNAL, FLXAU_INT, 'FLXAU', _RC)
      call MAPL_GetPointer(INTERNAL, FLCU_INT,  'FLCU',  _RC)
      call MAPL_GetPointer(INTERNAL, FLAU_INT,  'FLAU',  _RC)
      call MAPL_GetPointer(INTERNAL, FLXD_INT,  'FLXD',  _RC)
      call MAPL_GetPointer(INTERNAL, FLXAD_INT, 'FLXAD', _RC)
      call MAPL_GetPointer(INTERNAL, FLCD_INT,  'FLCD',  _RC)
      call MAPL_GetPointer(INTERNAL, FLAD_INT,  'FLAD',  _RC)
      call MAPL_GetPointer(INTERNAL, TS_INT,    'TS',    _RC)
      call MAPL_GetPointer(INTERNAL, DFDTS,     'DFDTS', _RC)
      call MAPL_GetPointer(INTERNAL, DFDTSC,    'DFDTSC',_RC)
      call MAPL_GetPointer(INTERNAL, DFDTSNA,   'DFDTSNA', _RC)
      call MAPL_GetPointer(INTERNAL, DFDTSCNA,  'DFDTSCNA',_RC)

      ! Determine calling sequence

      call MAPL_GetResource(MAPL,CalledLast,'CALLED_LAST:', default=1, _RC)

      ! Fill exported fluxed based on latest Ts

      if (CalledLast/=0) then
         call MAPL_TimerOn(MAPL,"-UPDATE_FLX")
         call Update_Flx(IM,JM,LM,_RC)
         call MAPL_TimerOff(MAPL,"-UPDATE_FLX")
      endif

      ! If it is time, refresh internal state.

      if (ESMF_AlarmIsRinging(ALARM,RC=STATUS)) then
         call ESMF_AlarmRingerOff(ALARM,_RC)

         call MAPL_TimerOn(MAPL,"-LW_DRIVER")
         call LW_Driver(IM,JM,LM,LATS,LONS,_RC)
         call MAPL_TimerOff(MAPL,"-LW_DRIVER")

      endif

      ! Fill exported fluxes based on latest Ts

      if (CalledLast==0) then
         call MAPL_TimerOn(MAPL,"-UPDATE_FLX")
         call Update_Flx(IM,JM,LM,_RC)
         call MAPL_TimerOff(MAPL,"-UPDATE_FLX")
      endif

      call MAPL_TimerOff(MAPL,"TOTAL")

      _RETURN(_SUCCESS)

   contains

      subroutine LW_Driver(IM,JM,LM,LATS,LONS,RC)

         ! RRTMGP module uses
         use mo_rte_kind,                only: wp
         use mo_gas_concentrations,      only: ty_gas_concs
         use mo_cloud_optics_rrtmgp,     only: ty_cloud_optics_rrtmgp
         use mo_cloud_sampling,          only: draw_samples, &
              sampled_mask_max_ran, sampled_mask_exp_ran, &
              sampled_urand_gen_max_ran
         use mo_optical_props,           only: ty_optical_props, &
              ty_optical_props_arry, ty_optical_props_1scl, &
              ty_optical_props_2str, ty_optical_props_nstr
         use mo_source_functions,        only: ty_source_func_lw
         use mo_fluxes,                  only: ty_fluxes_broadband
         use mo_fluxes_byband,           only: ty_fluxes_byband
         use mo_rte_lw,                  only: rte_lw
         use mo_load_coefficients,       only: load_and_init
         use mo_load_cloud_coefficients, only: load_cld_lutcoeff, load_cld_padecoeff

         ! used to avoid a reshaped copy of some arrays in RRTMGP blocking
         ! (Tom Clune's suggestion)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_loc, c_f_pointer

#ifdef HAVE_MKL
         ! Type of MKL VSL Basic RNGs
         ! (1) Mersenne Twister types
         ! brng = VSL_BRNG_MT19937
         ! Alternatives are VSL_BRNG_SFMT19937, maybe VSL_BRNG_MT2203?
         ! (2) Counter based PRNGs (CBPRNGs)
         ! brng = VSL_BRNG_PHILOX4X32X10  ! 10-round Philox 4x32 counter, 2x32 key
         ! Alternatives are VSL_BRNG_ARS5 ! faster if AES-NI instructions hardware supported
         !
         use MKL_VSL_TYPE
         use mo_rng_mklvsl_plus, only: ty_rng_mklvsl_plus
#else
         use mo_rng_mt19937, only: ty_rng_mt
#endif

         ! for RRTMGP (use implicit inside RRTMG)
         use cloud_condensate_inhomogeneity, only: condensate_inhomogeneous, zcw_lookup
         use cloud_subcol_gen, only : &
              correlation_length_cloud_fraction, correlation_length_condensate

         integer,                   intent(IN )    :: IM, JM, LM
         real,    dimension(IM,JM), intent(IN )    :: LATS, LONS
         integer, optional,         intent(OUT)    :: RC

         !  Locals

         character(len=ESMF_MAXSTR)        :: IAm
         integer                           :: STATUS
         integer                           :: loop_status

         ! local variables

         logical, parameter :: TRACE    = .true.

         integer, parameter :: NS       = 1       ! number of sub-grid surface types

         integer, parameter :: KICE     = 1
         integer, parameter :: KLIQUID  = 2
         integer, parameter :: KRAIN    = 3
         integer, parameter :: KSNOW    = 4
         integer, parameter :: KGRAUPEL = 5

         real    :: CO2_FIXED

         real    :: TAUCRIT                       ! pressure separating low and middle clouds
         real    :: PRS_LOW_MID                   ! pressure separating low and middle clouds
         real    :: PRS_MID_HIGH                  ! pressure separating low and high clouds
         integer :: LCLDMH                        ! model level separating high and middle clouds
         integer :: LCLDLM                        ! model level separating low  and middle clouds

         character(len=ESMF_MAXSTR), pointer :: AEROSOLS(:)

         integer :: i, j, K, L, YY, DOY, ibinary
         integer :: N !<<>> MSL

         real, dimension (IM,JM,NS)      :: FS    !  fractional cover of sub-grid regions
         real, dimension (IM,JM,NS)      :: TG    !  land or ocean surface temperature
         real, dimension (IM,JM,NS,10)   :: EG    !  land or ocean surface emissivity
         real, dimension (IM,JM,NS)      :: TV    !  vegetation temperature
         real, dimension (IM,JM,NS,10)   :: EV    !  vegetation emissivity
         real, dimension (IM,JM,NS,10)   :: RV    !  vegetation reflectivity
         real, dimension (IM,JM,LM,10)   :: TAUDIAG
         real, dimension (IM,JM,LM)      :: RH, PL, FCLD
         real, dimension (IM,JM,LM,5), target :: &
              CWC, &   ! in-cloud cloud water mixing ratio
              REFF     ! effective radius of cloud particles

         ! Local Aerosol Variables

         ! 4d: dimensioned (IM,JM,LM,NB_IRRAD)
         REAL, ALLOCATABLE, DIMENSION(:,:,:,:), target :: TAUA
         REAL, ALLOCATABLE, DIMENSION(:,:,:,:), target :: SSAA
         REAL, ALLOCATABLE, DIMENSION(:,:,:,:), target :: ASYA

         ! 3d pointer arrays to be associated with several 4d arrays for RRTMGP blocking.
         ! Collapses first two horizontal dimensions of these 4d arrays. OK since the
         ! latter arrays are not MAPL and so can be assumed contiguous.
         real, dimension(:,:,:), pointer :: TAUA_3d, SSAA_3d, ASYA_3d
         real, dimension(:,:,:), pointer :: CWC_3d, REFF_3d

         ! type(C_PTR) :: cptr  ! = c_loc(var), but done implicitly with c_loc below

         REAL :: X
         INTEGER :: IB, NA

         INTEGER :: OFFSET

         ! AERO state variables
         type (ESMF_State)                    :: AERO
         type (ESMF_Field)                    :: AS_FIELD
         character(len=ESMF_MAXSTR)           :: AS_FIELD_NAME
         type (ESMF_Field)                    :: AS_FIELD_Q
         integer                              :: AS_STATUS
         real, pointer,     dimension(:,:,:)  :: AS_PTR_3D
         real, pointer,     dimension(:,:,:)  :: AS_PTR_PLE
         real, pointer,     dimension(:,:,:)  :: AS_PTR_T
         real, pointer,     dimension(:,:,:)  :: AS_PTR_Q
         real, allocatable, dimension(:,:,:)  :: AS_ARR_RH
         real, allocatable, dimension(:,:,:)  :: AS_ARR_PL

         real, allocatable, dimension(:,:,:,:):: AEROSOL_EXT
         real, allocatable, dimension(:,:,:,:):: AEROSOL_SSA
         real, allocatable, dimension(:,:,:,:):: AEROSOL_ASY

         real, pointer,     dimension(:,:,:)  :: VAR_PTR_3D

         logical                              :: implements_aerosol_optics

         integer                              :: band

         ! Variables for RRTMG Code

         integer :: iceflglw        ! Flag for ice particle specification
         integer :: liqflglw        ! Flag for liquid droplet specification
         logical :: Ts_derivs       ! calculate Tsurf derivatives of upward fluxes
         integer :: NN, IJ, LV

         real,    allocatable, dimension(:,:)   :: FCLD_R
         real,    allocatable, dimension(:,:)   :: TLEV_R       ! Edge Level temperature
         real,    allocatable, dimension(:,:)   :: PLE_R        ! Reverse of level pressure
         real,    allocatable, dimension(:,:)   :: ZM_R         ! Reverse of layer height
         real,    allocatable, dimension(:,:)   :: EMISS        ! Surface emissivity at 16 RRTMG bands
         real,    allocatable, dimension(:,:)   :: CLIQWP       ! Cloud liquid water path
         real,    allocatable, dimension(:,:)   :: CICEWP       ! Cloud ice water path
         real,    allocatable, dimension(:,:)   :: RELIQ        ! Cloud liquid effective radius
         real,    allocatable, dimension(:,:)   :: REICE        ! Cloud ice effective radius
         real,    allocatable, dimension(:,:,:) :: TAUAER
         real,    allocatable, dimension(:,:)   :: PL_R, T_R,  Q_R, O2_R,  O3_R
         real,    allocatable, dimension(:,:)   :: CO2_R, CH4_R, N2O_R, CFC11_R, CFC12_R, CFC22_R, CCL4_R
         real,    allocatable, dimension(:)     :: TSFC
         real,    allocatable, dimension(:,:)   :: UFLX, DFLX, UFLXC, DFLXC, DUFLX_DTS, DUFLXC_DTS
         integer, allocatable, dimension(:,:)   :: CLEARCOUNTS
         real,    allocatable, dimension(:)     :: ALAT
         real,    allocatable, dimension(:,:)   :: OLRBRG, DOLRBRG_DTS

         ! pmn: should we update these?
         real, parameter :: O2   = 0.2090029E+00 ! preexisting
         real, parameter :: N2   = 0.7906400E+00 ! approx from rrtmgp input file
         real, parameter :: CCL4 = 0.1105000E-09 ! preexisting
         real, parameter :: CO   = 0.            ! currently zero

         ! variables for RRTMGP code

         ! conversion factor (see below)
         real(wp), parameter :: cwp_fac = real(1000./MAPL_GRAV,kind=wp)

         ! input arrays: dimensions (ncol, nlay[+1]) [Pa,K]
         real(wp), dimension(:,:), allocatable         :: p_lay, t_lay, dp_wp, cf_wp
         real(wp), dimension(:,:), allocatable         :: p_lev
         real(wp), dimension(:,:), allocatable, target :: t_lev

         ! inter-layer separations (from mid-points) (ncol,nlay-1) [m]
         real(wp), dimension(:,:), allocatable         :: dzmid

         ! surface input arrays
         real(wp), dimension(:),   allocatable         :: t_sfc
         real(wp), dimension(:,:), allocatable         :: emis_sfc ! first dim is band

         ! fluxes:
         ! broadband
         real(wp), dimension(:,:), allocatable, target :: &
              flux_up_clrsky, flux_dn_clrsky, dfupdts_clrsky, &
              flux_up_clrnoa, flux_dn_clrnoa, dfupdts_clrnoa, &
              flux_up_allsky, flux_dn_allsky, dfupdts_allsky, &
              flux_up_allnoa, flux_dn_allnoa, dfupdts_allnoa
         ! byband
         real(wp), dimension(:,:,:), allocatable, target :: &
              bnd_flux_up_allnoa, bnd_dfupdts_allnoa, &
              bnd_flux_up_allsky, bnd_dfupdts_allsky

         ! derived types for interacting with RRTMGP
         type(ty_gas_optics_rrtmgp), pointer           :: k_dist
         type(ty_gas_concs)                            :: gas_concs, gas_concs_block
         type(ty_cloud_optics_rrtmgp)                  :: cloud_optics
         type(ty_source_func_lw)                       :: sources
         type(ty_fluxes_broadband)                     :: fluxes_clrsky, fluxes_clrnoa, fluxes_allnoa, fluxes_allsky
         type(ty_fluxes_byband)                        :: fluxes_byband_allnoa, fluxes_byband_allsky

         ! The band-space (ncols_block,nlay,nbnd) aerosol and in-cloud optical properties
         ! Polymorphic with dynamic type (#streams) defined later
         class(ty_optical_props_arry), allocatable :: cloud_props_bnd, aer_props

         ! The g-point cloud optical properties used for mcICA
         class(ty_optical_props_arry), allocatable :: cloud_props_gpt

         ! The g-point optical properties used in RT calculations for clean|dirty exports
         ! Polymorphic with dynamic type (#streams) defined later
         class(ty_optical_props_arry), allocatable :: clean_optical_props, dirty_optical_props

         ! RRTMGP locals
         logical :: top_at_1, u2s, partial_block, gen_mro, cond_inhomo
         logical :: need_dirty_optical_props, need_cloud_optical_props
         logical :: export_clrnoa, export_clrsky, export_allnoa, export_allsky
         logical ::   calc_clrnoa,   calc_clrsky,   calc_allnoa,   calc_allsky
         logical :: allnoa_to_allsky_band_xfer_needed
         integer :: ncol, nbnd, ngpt, nmom, nga, icergh
         integer :: b, nBlocks, colS, colE, ncols_block, &
              partial_blockSize, icol, isub, ilay, igpt
         character(len=ESMF_MAXPATHLEN) :: k_dist_file, cloud_optics_file
         character(len=ESMF_MAXSTR)     :: error_msg
         character(len=128)             :: cloud_optics_type, cloud_overlap_type
         type (ESMF_Time)               :: ReferenceTime
         type (ESMF_TimeInterval)       :: RefreshInterval
         real(wp) :: cld_frac
         real :: sigma_qcw

         ! for global gcolumn index seeding of PRNGs
         integer :: iBeg, iEnd, jBeg, jEnd
         integer :: IM_World, JM_World, Gdims(3)

         ! a column random number generator
#ifdef HAVE_MKL
         type(ty_rng_mklvsl_plus) :: rng
#else
         type(ty_rng_mt) :: rng
#endif
         integer, dimension(:), allocatable :: seeds

         ! uniform random numbers need by mcICA (ngpt,nlay,rrtmgp_blocksize)
         real(wp), dimension(:,:,:), allocatable :: &
              urand, urand_aux, urand_cond, urand_cond_aux

         ! Cloud mask for overlap scheme (ncols_block,nlay,ngpt)
         logical,  dimension(:,:,:), allocatable :: cld_mask

         ! sub-gridscale condensate scaling for overlap scheme (ncols_block,nlay,ngpt)
         real(wp), dimension(:,:,:), allocatable :: zcw

         ! correlation length scales [m] for cloud presence and condensate (ncol)
         real, dimension(:), allocatable :: adl, rdl

         ! binomial probability of maximum overlap (cf. random overlap)
         ! for cloud presence and condensate (ncols_block,nlay-1)
         real(wp), dimension(:,:), allocatable :: alpha, rcorr

         ! TEMP ... see below
         real(wp) :: press_ref_min, ptop
         real(wp) ::  temp_ref_min, tmin
         real(wp) ::  temp_ref_max, tmax

         ! block size for efficient column processing (set from resource file)
         integer :: rrtmgp_blockSize

         ! For aerosol
         integer                    :: in
         real                       :: xx, LWT, IWT
         type (ESMF_Time)           :: CURRENTTIME
         real, dimension (LM+1)     :: TLEV
         real, dimension (LM)       :: DP

         ! pointers to import

         real, pointer, dimension(:    )   :: PREF
         real, pointer, dimension(:,:  )   :: TS
         real, pointer, dimension(:,:  )   :: EMIS
         real, pointer, dimension(:,:,:)   :: PLE, T,  Q,  O3
         real, pointer, dimension(:,:,:)   :: CH4, N2O, CFC11, CFC12, HCFC22
         real, pointer, dimension(:,:,:)   :: QL, QI, QR, QS, QG
         real, pointer, dimension(:,:,:)   :: RI, RL, RR, RS, RG, FCLD_IN
         real, pointer, dimension(:,:,:,:) :: RAERO
         real, pointer, dimension(:,:,:)   :: QAERO
         real, pointer, dimension(:,:,:)   :: CO2_3d => null() ! <<>> MSL
         real, pointer, dimension(:,:,:)   :: tmp_3d => null() ! <<>> MSL

         ! pointers to exports

         real, pointer, dimension(:,:  )   :: CLDPRS
         real, pointer, dimension(:,:  )   :: CLDTMP
         real, pointer, dimension(:,:,:)   :: TAUIR
         real, pointer, dimension(:,:  )   :: CLDTTLW
         real, pointer, dimension(:,:  )   :: CLDHILW
         real, pointer, dimension(:,:  )   :: CLDMDLW
         real, pointer, dimension(:,:  )   :: CLDLOLW
         real, pointer, dimension(:,:  )   :: TSREFF
         real, pointer, dimension(:,:  )   :: SFCEM
         real, pointer, dimension(:,:  )   :: LWS0
         real, pointer, dimension(:,:  )   :: DSFDTS

         ! for compact multi-export handling
         real, pointer, dimension(:,:  ) :: ptr2d
         real, pointer, dimension(:,:,:) :: ptr3d

         type(StringVector) :: string_vec
         type(StringVectorIterator) :: string_vec_iter
         character(len=:), pointer :: string_pointer

         ! helper for testing RRTMGP error status on return;
         ! allows line number reporting cf. original call method
#define TEST_(A) error_msg = A; if (trim(error_msg)/="") then; _ASSERT(.false.,"RRTMGP Error: "//trim(error_msg)); endif

         logical :: USE_PRECIP_IN_RADIATION
         integer :: PARTITION_SIZE

         real, parameter :: SSA_MAX = 0.999999
         real, parameter :: ASY_MAX = 0.999

         !  Begin...

         IAm = "LW_Driver"
         call MAPL_TimerOn(MAPL,"--MISC")

         ! Pointer to Imports used only for full transfer calculation

         call MAPL_GetPointer(IMPORT, PLE,    'PLE',    _RC)
         call MAPL_GetPointer(IMPORT, T,      'T',      _RC)
         call MAPL_GetPointer(IMPORT, Q,      'QV',     _RC)
         call MAPL_GetPointer(IMPORT, QL,     'QL',     _RC)
         call MAPL_GetPointer(IMPORT, QI,     'QI',     _RC)
         call MAPL_GetPointer(IMPORT, QR,     'QR',     _RC)
         call MAPL_GetPointer(IMPORT, QS,     'QS',     _RC)
         call MAPL_GetPointer(IMPORT, QG,     'QG',     _RC)
         call MAPL_GetPointer(IMPORT, RL,     'RL',     _RC)
         call MAPL_GetPointer(IMPORT, RI,     'RI',     _RC)
         call MAPL_GetPointer(IMPORT, RR,     'RR',     _RC)
         call MAPL_GetPointer(IMPORT, RS,     'RS',     _RC)
         call MAPL_GetPointer(IMPORT, RG,     'RG',     _RC)
         call MAPL_GetPointer(IMPORT, O3,     'O3',     _RC)
         call MAPL_GetPointer(IMPORT, CH4,    'CH4',    _RC)
         call MAPL_GetPointer(IMPORT, N2O,    'N2O',    _RC)
         call MAPL_GetPointer(IMPORT, CFC11,  'CFC11',  _RC)
         call MAPL_GetPointer(IMPORT, CFC12,  'CFC12',  _RC)
         call MAPL_GetPointer(IMPORT, HCFC22, 'HCFC22', _RC)
         call MAPL_GetPointer(IMPORT, FCLD_IN,'FCLD',   _RC)
         call MAPL_GetPointer(IMPORT, EMIS,   'EMIS',   _RC)
         call MAPL_GetPointer(IMPORT, PREF,   'PREF',   _RC)
         call MAPL_GetPointer(IMPORT, TS,     'TS',     _RC)

         PL = 0.5*(PLE(:,:,:UBOUND(PLE,3)-1)+PLE(:,:,LBOUND(PLE,3)+1:))
         RH = Q/GEOS_QSAT(T,PL,PASCALS=.true.)

         ! make a copy of 'FCLD' so can optionally change it without changing import state
         FCLD = FCLD_IN

         ! Option to force binary clouds for LW
         call MAPL_GetResource(MAPL,ibinary,"RADLW_BINARY_CLOUDS:",DEFAULT=0,_RC)
         if (ibinary /= 0) where (FCLD > 0.) FCLD = 1.

         ! Get trace gases concentrations by volume (pppv) from configuration

         call MAPL_GetResource (MAPL, CO2_FIXED, 'CO2:', _RC)

         ! <<>> MSL
         if(CO2_FIXED.eq.-2.0) then ! 3D CO2
            if (USE_CHOU) then ! No 3D CO2 if USE_CHOU
               CO2_FIXED = -1.0
            else
               call MAPL_GetPointer(IMPORT, CO2_3d, 'CO2', _RC)
               call ESMF_ClockGet(CLOCK, currTIME=CURRENTTIME, _RC)
               call ESMF_TimeGet (CURRENTTIME, YY=YY, DayOfYear=DOY, _RC)
               CO2_FIXED = GETCO2(YY,DOY)
               call MAPL_GetPointer(EXPORT, tmp_3d, 'CO2_FIXED', NotFoundOK=.true., RC=STATUS)
               if (associated(tmp_3d)) then
                  tmp_3d = CO2_FIXED
                  tmp_3d => null()
               endif
            endif
         endif

         if(CO2_FIXED.eq.-1.0) then
            call ESMF_ClockGet(CLOCK, currTIME=CURRENTTIME, _RC)
            call ESMF_TimeGet (CURRENTTIME, YY=YY, DayOfYear=DOY, _RC)
            CO2_FIXED = GETCO2(YY,DOY)
         endif

         call MAPL_GetResource (MAPL, PRS_LOW_MID,  'PRS_LOW_MID_CLOUDS:',  DEFAULT=70000., _RC)
         call MAPL_GetResource (MAPL, PRS_MID_HIGH, 'PRS_MID_HIGH_CLOUDS:', DEFAULT=40000., _RC)

         ! Set up the RATS toggles  <<>> MSL
         ! -- these fields will be turn on/off to eval flux impacts
         ! -- ideally, we could query the exports to find if any actually -need- computing
         !    because if not (e.g. CO2 is listed as a RAT_DIAG, but HISTORY.rc has
         !    no diagnostic output for that RAT), there's no need to run an additional RRTMG_LW().
         ! -- This is done every call to Run(), when it really only needs to be done once
         IF (first) then

            call ESMF_ConfigFindLabel(CF, 'RATS_DIAGNOSTICS:', RC=STATUS) ! Use STATUS to test if label was found

            nRATS = 0 ! Default, no RAT diags

            ! No error thrown. Just go around this if nothing learnable from config.
            IF (STATUS .eq. ESMF_SUCCESS) THEN ! if the label was found...

               ! Get number of words in config line
               nRATS = ESMF_ConfigGetLen(CF,label='RATS_DIAGNOSTICS:',_RC)

               allocate(nameRATS(nRATS), STAT=STATUS)
               VERIFY_(STATUS)

               ! Put the cursor at the label
               call ESMF_ConfigFindLabel(CF, 'RATS_DIAGNOSTICS:', _RC)

               DO i=1,nRATS
                  call ESMF_ConfigGetAttribute(CF,gen_str,_RC)
                  nameRATS(i) = trim(gen_str)
               ENDDO

               ! Only allocate this if needed
               allocate(TMP_R(IM*JM,LM),_STAT)
            ENDIF
            first = .false. ! Don't repeat this.
         ENDIF ! first

         ! Prepare for aerosol optics calculations

         ! Set the offset for the IRRAD aerosol bands
         if (USE_RRTMGP_SORAD) then
            OFFSET = NB_RRTMGP_SORAD
         else if (USE_RRTMG_SORAD) then
            OFFSET = NB_RRTMG_SORAD
         else
            OFFSET = NB_CHOU_SORAD
         end if

         ! For now, use the same emissivity for all bands

         do K = 1, 10
            EG(:,:,1,K)   = EMIS(:,:)
         end do

         ! For now, hardwire vegetation and aerosol parameters

         FS                  = 1.0
         TG(:,:,1)           = TS
         TV(:,:,1)           = TS
         EV                  = 0.0
         RV                  = 0.0

         ! Copy cloud constituent properties into contiguous buffers

         ! In-cloud water contents
         CWC (:,:,:,KICE    ) = QI
         CWC (:,:,:,KLIQUID ) = QL
         CWC (:,:,:,KRAIN   ) = QR
         CWC (:,:,:,KSNOW   ) = QS
         CWC (:,:,:,KGRAUPEL) = QG

         ! Effective radii [microns]
         REFF(:,:,:,KICE    ) = RI * 1.0e6
         REFF(:,:,:,KLIQUID ) = RL * 1.0e6
         REFF(:,:,:,KRAIN   ) = RR * 1.0e6
         REFF(:,:,:,KSNOW   ) = RS * 1.0e6
         REFF(:,:,:,KGRAUPEL) = RG * 1.0e6
         WHERE (RI == MAPL_UNDEF) REFF(:,:,:,KICE    ) = 36.
         WHERE (RL == MAPL_UNDEF) REFF(:,:,:,KLIQUID ) = 14.
         WHERE (RR == MAPL_UNDEF) REFF(:,:,:,KRAIN   ) = 50.
         WHERE (RS == MAPL_UNDEF) REFF(:,:,:,KSNOW   ) = 50.
         WHERE (RG == MAPL_UNDEF) REFF(:,:,:,KGRAUPEL) = 50.

         ! Determine the model level separating high-middle and low-middle clouds

         _ASSERT(PRS_MID_HIGH > PREF(1)     , 'mid-high pressure band boundary too high!')
         _ASSERT(PRS_LOW_MID  > PRS_MID_HIGH, 'pressure band misordering!')
         _ASSERT(PRS_LOW_MID  < PREF(LM)    , 'low-mid pressure band boundary too low!')

         ! find mid-high interface level
         k = 1
         do while ( PREF(k) < PRS_MID_HIGH )
            k=k+1
         end do
         LCLDMH = k
         ! Guaranteed that LCLDMH > 1 (by first ASSERT above)
         !    and that PREF(LCLDMH) >= PRS_MID_HIGH (by while loop)

         ! find low-mid interface level
         do while ( PREF(k) < PRS_LOW_MID )
            k=k+1
         end do
         LCLDLM = k
         ! Guaranteed that LCLDLM <= LM (by third assert above)
         !    and that PREF(LCLDLM) >= PRS_LOW_MID (by while loop)

         ! But it's still possible that LCLDLM == LCLDMH if the
         ! interface pressures are too close. We now ASSERT to
         ! prevent this.
         _ASSERT(LCLDMH < LCLDLM, 'PRS_LOW_MID and PRS_MID_HIGH are too close!')

         ! now we have 1 < LCLDMH < LCLDLM <= LM and can use:
         !    layers [1,      LCLDMH-1] are in high pressure band
         !    layers [LCLDMH, LCLDLM-1] are in mid  pressure band
         !    layers [LCLDLM, LM      ] are in low  pressure band

         call MAPL_GetPointer(EXPORT, CLDTTLW, 'CLDTTLW', _RC)
         call MAPL_GetPointer(EXPORT, CLDHILW, 'CLDHILW', _RC)
         call MAPL_GetPointer(EXPORT, CLDMDLW, 'CLDMDLW', _RC)
         call MAPL_GetPointer(EXPORT, CLDLOLW, 'CLDLOLW', _RC)

         ! Begin aerosol code

         ! Allocate per-band aerosol arrays

         ALLOCATE (TAUA(IM,JM,LM,NB_IRRAD),_STAT)
         ALLOCATE (SSAA(IM,JM,LM,NB_IRRAD),_STAT)
         ALLOCATE (ASYA(IM,JM,LM,NB_IRRAD),_STAT)

         ! Zero out aerosol arrays. If NA == 0, these zeroes are then used inside IRRAD.
         NA   = 0

         TAUA = 0.
         SSAA = 0.
         ASYA = 0.

         ! If we have aerosols, accumulate the arrays

         call MAPL_TimerOn(MAPL,"---AEROSOLS")

         call ESMF_StateGet(IMPORT, 'AERO', AERO, _RC)

         call ESMF_AttributeGet(aero, name='implements_aerosol_optics_method', &
              value=implements_aerosol_optics, _RC)

         RADIATIVELY_ACTIVE_AEROSOLS: if (implements_aerosol_optics) then

            ! set RH for aerosol optics
            call ESMF_AttributeGet(AERO, name='relative_humidity_for_aerosol_optics', value=AS_FIELD_NAME, _RC)

            if (AS_FIELD_NAME /= '') then
               call MAPL_GetPointer(AERO, AS_PTR_3D, trim(AS_FIELD_NAME), _RC)

               AS_PTR_3D = RH
            end if

            ! set PLE for aerosol optics
            call ESMF_AttributeGet(AERO, name='air_pressure_for_aerosol_optics', value=AS_FIELD_NAME, _RC)

            if (AS_FIELD_NAME /= '') then
               call MAPL_GetPointer(AERO, AS_PTR_3D, trim(AS_FIELD_NAME), _RC)

               AS_PTR_3D = PLE
            end if

            ! allocate memory for total aerosol ext, ssa and asy at all solar bands
            allocate(AEROSOL_EXT(IM,JM,LM,NB_IRRAD),  &
                 AEROSOL_SSA(IM,JM,LM,NB_IRRAD),  &
                 AEROSOL_ASY(IM,JM,LM,NB_IRRAD),  stat=STATUS)
            VERIFY_(STATUS)

            AEROSOL_EXT = 0.
            AEROSOL_SSA = 0.
            AEROSOL_ASY = 0.

            ! compute aerosol optics at all solar bands
            IR_BANDS: do band = 1, NB_IRRAD
               call ESMF_AttributeSet(AERO, name='band_for_aerosol_optics', value=(OFFSET+band), _RC)

               ! execute the aero provider's optics method
               call ESMF_MethodExecute(AERO, label="run_aerosol_optics", _RC)

               ! EXT from AERO_PROVIDER
               call ESMF_AttributeGet(AERO, name='extinction_in_air_due_to_ambient_aerosol', value=AS_FIELD_NAME, _RC)

               if (AS_FIELD_NAME /= '') then
                  call MAPL_GetPointer(AERO, AS_PTR_3D, trim(AS_FIELD_NAME),  _RC)

                  if (associated(AS_PTR_3D)) then
                     AEROSOL_EXT(:,:,:,band) = MAX(AS_PTR_3D,0.0)
                  end if
               end if

               ! SSA from AERO_PROVIDER
               call ESMF_AttributeGet(AERO, name='single_scattering_albedo_of_ambient_aerosol', value=AS_FIELD_NAME, _RC)

               if (AS_FIELD_NAME /= '') then
                  call MAPL_GetPointer(AERO, AS_PTR_3D, trim(AS_FIELD_NAME),  _RC)

                  if (associated(AS_PTR_3D)) then
                     AEROSOL_SSA(:,:,:,band) = MIN(MAX(AS_PTR_3D,0.0),SSA_MAX)
                  end if
               end if

               ! ASY from AERO_PROVIDER
               call ESMF_AttributeGet(AERO, name='asymmetry_parameter_of_ambient_aerosol', value=AS_FIELD_NAME, _RC)

               if (AS_FIELD_NAME /= '') then
                  call MAPL_GetPointer(AERO, AS_PTR_3D, trim(AS_FIELD_NAME),  _RC)

                  if (associated(AS_PTR_3D)) then
                     AEROSOL_ASY(:,:,:,band) = MIN(MAX(AS_PTR_3D,0.0),ASY_MAX)
                  end if
               end if
            end do IR_BANDS

            NA = 3

            TAUA = AEROSOL_EXT
            SSAA = AEROSOL_SSA
            ASYA = AEROSOL_ASY

            deallocate(AEROSOL_EXT, _STAT)
            deallocate(AEROSOL_SSA, _STAT)
            deallocate(AEROSOL_ASY, _STAT)

         end if RADIATIVELY_ACTIVE_AEROSOLS

         call MAPL_TimerOff(MAPL,"---AEROSOLS")

         call MAPL_TimerOff(MAPL,"--MISC")

         SCHEME: if (USE_CHOU) then

            call MAPL_TimerOn (MAPL,"--IRRAD",_RC)

            ! Do longwave calculations on a list of soundings
            !  This fills the internal state
            ! Note: IRRAD wants all species in mole fraction
            ! except O3, which must be in mass mixing ratio.

            call MAPL_TimerOn(MAPL,"---IRRAD_RUN",_RC)
            call IRRAD( IM*JM, LM,       PLE,                           &
                 T,        Q,      O3,    TS,     CO2_FIXED,                &
                 TRACE,    N2O,   CH4,    CFC11,     CFC12, HCFC22,         &
                 CWC,    FCLD,  LCLDMH, LCLDLM,    REFF,                    &
                 NS,       FS,     TG,    EG,     TV,        EV,    RV,     &
                 NA, NB_CHOU, TAUA, SSAA, ASYA,                             &
                 FLXU_INT,  FLCU_INT, FLAU_INT, FLXAU_INT,                  &
                 FLXD_INT,  FLCD_INT, FLAD_INT, FLXAD_INT,                  &
                 DFDTS, SFCEM_INT, TAUDIAG                                  )
            call MAPL_TimerOff(MAPL,"---IRRAD_RUN",_RC)

            ! pmn:
            ! Chou-Suarez does not provide these derivatives
            ! so clear is set to zero, no-aerosol to aerosol
            DFDTSC = 0.
            DFDTSNA  = DFDTS
            DFDTSCNA = DFDTSC

            call MAPL_TimerOff(MAPL,"--IRRAD",_RC)

         else if (USE_RRTMGP) then

            call MAPL_TimerOn(MAPL,"--RRTMGP",_RC)

            ! columns are independent so collapse horizontal to 1D
            ncol = IM*JM

            ! absorbing gas names
            error_msg = gas_concs%init([character(3) :: &
                 'h2o','co2','o3','n2o','co','ch4','o2','n2'])
            TEST_(error_msg)

            if (associated(  CO2_3d)) &
                 allocate(CO2_R(IM*JM,LM),_STAT)
            allocate(  Q_R(IM*JM,LM),_STAT)
            allocate( O3_R(IM*JM,LM),_STAT)
            allocate(N2O_R(IM*JM,LM),_STAT)
            allocate(CH4_R(IM*JM,LM),_STAT)

            if (associated(  CO2_3d)) &
                 CO2_R = reshape( CO2_3d                          ,(/ncol,LM/))
            Q_R = reshape( Q/(1.-Q)*(MAPL_AIRMW/MAPL_H2OMW),(/ncol,LM/))
            O3_R = reshape( O3      *(MAPL_AIRMW/MAPL_O3MW ),(/ncol,LM/))
            N2O_R = reshape( N2O                             ,(/ncol,LM/))
            CH4_R = reshape( CH4                             ,(/ncol,LM/))

            ! Clean up negatives
            if (associated(  CO2_3d)) &
                 WHERE ( CO2_R < 0.) CO2_R = 0.
            WHERE (   Q_R < 0.)   Q_R = 0.
            WHERE (  O3_R < 0.)  O3_R = 0.
            WHERE ( N2O_R < 0.) N2O_R = 0.
            WHERE ( CH4_R < 0.) CH4_R = 0.

            ! load gas concentrations (volume mixing ratios)
            ! "constant" gases
            TEST_(gas_concs%set_vmr('n2' , real(N2 ,kind=wp)))
            TEST_(gas_concs%set_vmr('o2' , real(O2 ,kind=wp)))
            if (.not. associated(CO2_3d)) TEST_(gas_concs%set_vmr('co2', real(CO2_FIXED,kind=wp))) ! <<>> MSL
            TEST_(gas_concs%set_vmr('co' , real(CO ,kind=wp)))
            ! variable gases
            ! (ozone converted from mass mixing ratio, water vapor from specific humidity)
            if (associated(  CO2_3d)) then
               TEST_(gas_concs%set_vmr('co2', real(CO2_R,kind=wp)))
            else
               TEST_(gas_concs%set_vmr('co2', real(CO2_FIXED,kind=wp))) ! <<>> MSL
            endif
            TEST_(gas_concs%set_vmr('h2o', real(  Q_R,kind=wp)))
            TEST_(gas_concs%set_vmr('o3' , real( O3_R,kind=wp)))
            TEST_(gas_concs%set_vmr('n2o', real(N2O_R,kind=wp)))
            TEST_(gas_concs%set_vmr('ch4', real(CH4_R,kind=wp)))
            if (associated(CO2_3d)) TEST_(gas_concs%set_vmr('co2', real(reshape(CO2_3d  ,(/ncol,LM/)),kind=wp))) !<<>> MSL

            if (associated(  CO2_3d)) &
                 deallocate( CO2_R,_STAT)
            deallocate(   Q_R,_STAT)
            deallocate(  O3_R,_STAT)
            deallocate( N2O_R,_STAT)
            deallocate( CH4_R,_STAT)

            ! access RRTMGP internal state from the GC
            call ESMF_UserCompGetInternalState(GC, 'RRTMGP_state', wrap, status)
            VERIFY_(status)
            rrtmgp_state => wrap%ptr

            ! initialize k-distribution if not already done
            call MAPL_GetResource( &
                 MAPL, k_dist_file, "RRTMGP_GAS_LW:", &
                 DEFAULT='rrtmgp-gas-lw-g128.nc',_RC)
            if (.not. rrtmgp_state%initialized) then
               ! gas_concs needed only to access required gas names
               call MAPL_TimerOn(MAPL,"---RRTMGP_IO_GAS",_RC)
               call load_and_init(rrtmgp_state%k_dist, trim(k_dist_file), gas_concs)
               call MAPL_TimerOff(MAPL,"---RRTMGP_IO_GAS",_RC)
               if (.not. rrtmgp_state%k_dist%source_is_internal()) then
                  TEST_("RRTMGP-LW: does not seem to be LW")
               endif
               rrtmgp_state%initialized = .true.
            endif

            ! access by shorter name
            k_dist => rrtmgp_state%k_dist

            ! spectral dimensions
            ngpt = k_dist%get_ngpt()
            nbnd = k_dist%get_nband()
            _ASSERT(nbnd == NB_RRTMGP, 'RRTMGP-LW: expected different number of bands')

            ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! For reference, comparison of RRTMG and RRTMGP bands:
            ! from RRTMG:
            ! wavenum1(:) = (/ 10., 350., 500., 630., 700., 820.,  980., 1080., 1180., 1390., 1480., 1800., 2080., 2250., 2380., 2600./)
            ! wavenum2(:) = (/350., 500., 630., 700., 820., 980., 1080., 1180., 1390., 1480., 1800., 2080., 2250., 2380., 2600., 3250./)
            ! from RRTMGP:
            ! write(*,*) 'band_lims_wvn(2,nbnd):', k_dist%get_band_lims_wavenumber()  ! with output reordered
            !                  10., 250., 500., 630., 700., 820.,  980., 1080., 1180., 1390., 1480., 1800., 2080., 2250., 2390., 2680.
            !                 250., 500., 630., 700., 820., 980., 1080., 1180., 1390., 1480., 1800., 2080., 2250., 2390., 2680., 3250.
            ! clearly there are some differences (250, 2390, 2680) ... have redone aerosol tables
            ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            ! allocate input arrays
            allocate(t_sfc(ncol), emis_sfc(nbnd,ncol), _STAT)
            allocate(p_lay(ncol,LM), t_lay(ncol,LM), dp_wp(ncol,LM), cf_wp(ncol,LM), _STAT)
            allocate(p_lev(ncol,LM+1), t_lev(ncol,LM+1), _STAT)
            allocate(dzmid(ncol,LM-1), _STAT)

            ! load input arrays ...
            ! surface properties
            ! for now, use the same emissivity for all bands
            t_sfc    = real(       reshape(TS  ,(/ncol/))        ,kind=wp)
            emis_sfc = real(spread(reshape(EMIS,(/ncol/)),1,nbnd),kind=wp)

            ! Currently k_dist%temp_ref_max = 355K ~ 82C, but GEOS-5 seems to
            ! sometimes exceed the maximum temperature. See more comments under
            ! layer temperature kluge below. We clip it here as a kluge.
            temp_ref_max = k_dist%get_temp_max() - 0.01_wp
            tmax = maxval(t_sfc)
            where (t_sfc > temp_ref_max) t_sfc = temp_ref_max

            ! basic profiles
            p_lay = real(reshape(PL  ,(/ncol,LM  /)), kind=wp)
            t_lay = real(reshape(T   ,(/ncol,LM  /)), kind=wp)
            p_lev = real(reshape(PLE ,(/ncol,LM+1/)), kind=wp)
            cf_wp = real(reshape(FCLD,(/ncol,LM  /)), kind=wp)

            ! RRTMGP's rte_lw takes a vertical ordering flag
            ! (no need to flip columns as with RRTMG)
            top_at_1 = p_lay(1, 1) < p_lay(1, LM)
            _ASSERT(top_at_1, 'unexpected vertical ordering')

            ! layer pressure thicknesses used for cloud water path calculations
            ! (do before any KLUGE to top pressure so optical paths wont be affected)
            ! (also better to use these unKLUGED pressure intervals in t_lev calculation)
            dp_wp = p_lev(:,2:LM+1) - p_lev(:,1:LM)

            ! Because currently k_dist%press_ref_min ~ 1.005 > GEOS-5 ptop of 1.0 Pa.
            ! Find better solution, perhaps getting AER to add a higher top.
            press_ref_min = k_dist%get_press_min()
            where (p_lev(:,1) < press_ref_min) p_lev(:,1) = press_ref_min
            ! make sure no pressure ordering issues were created
            _ASSERT(all(p_lev(:,1) < p_lay(:,1)), 'pressure kluge causes misordering')

            ! pmn: temperature KLUGE
            ! Find better solution, perhaps getting AER to produce a table with a
            ! larger temperature range.
            temp_ref_min = k_dist%get_temp_min() + 0.01_wp
            where (t_lay < temp_ref_min) t_lay = temp_ref_min
            temp_ref_max = k_dist%get_temp_max() - 0.01_wp
            where (t_lay > temp_ref_max) t_lay = temp_ref_max

            ! Calculate interface temperatures (t_lev) and layer midpoint separations (dzmid)
            ! pmn: t_lev is an optional argument of gas_optics(), and if not provided, it will supply its
            !   own internally. Could try running with this latter option to see what difference it makes.
            ! pmn: these t_lev must also be >= temp_ref_min. Since the core of the t_lev calculation below
            !   is an INTERPOLATION, and since the t_lay are already KLUGED to >= temp_ref_min, this should
            !   not be a problem. But this is why the t_lev calculation must occur AFTER the t_lay KLUGE.
            !   Note that t_lev(1) gets a copy of t_lev(2), so will also be in range. We are not worried
            !   about TS being < temp_ref_min = 160K (surface values wont get that cold!)
            ! dzmid(k) is separation [m] between midpoints of layers k and k+1 (sign not important, positive
            !   here). dz ~ RT/g x dp/p by hydrostatic eqn and ideal gas eqn. The jump from LAYER k to k+1
            !   is centered on LEVEL k+1 since the LEVEL indices are one-based.
            do k = 1,LM-1
               ! t_lev interpolated between neighboring t_lay
               t_lev(:,k+1) = (t_lay(:,k) * dp_wp(:,k+1) + t_lay(:,k+1) * dp_wp(:,k)) / (dp_wp(:,k+1) + dp_wp(:,k))
               dzmid(:,k) = t_lev(:,k+1) * real(MAPL_RGAS/MAPL_GRAV,kind=wp) * (p_lay(:,k+1) - p_lay(:,k)) / p_lev(:,k+1)
            end do
            t_lev(:,1) = t_lev(:,2)                              ! assume isotropic at TOA
            t_lev(:,LM+1) = real(reshape(TS,(/ncol/)),kind=wp)  ! ~surface air temperature

            ! for efficiency sake, we try to calculate only what we export ...

            ! are clear clean exports requested?
            export_clrnoa = .false.

            call string_vec%push_back('FLA')
            call string_vec%push_back('FLAD')
            call string_vec%push_back('FLAU')
            string_vec_iter = string_vec%begin()
            do while ( string_vec_iter /= string_vec%end() )
               string_pointer => string_vec_iter%get()
               call MAPL_GetPointer( EXPORT, ptr3d, string_pointer, _RC)
               export_clrnoa = (export_clrnoa .or. associated(ptr3d))
               call string_vec_iter%next()
            end do

            call string_vec%clear()
            call string_vec%push_back('OLA')
            call string_vec%push_back('FLNSA')
            call string_vec%push_back('LAS')
            string_vec_iter = string_vec%begin()
            do while ( string_vec_iter /= string_vec%end() )
               string_pointer => string_vec_iter%get()
               call MAPL_GetPointer( EXPORT, ptr2d, string_pointer, _RC)
               export_clrnoa = (export_clrnoa .or. associated(ptr2d))
               call string_vec_iter%next()
            end do

            ! are clear dirty exports requested?
            export_clrsky = .false.

            call string_vec%clear()
            call string_vec%push_back('FLC')
            call string_vec%push_back('FLCD')
            call string_vec%push_back('FLCU')
            string_vec_iter = string_vec%begin()
            do while ( string_vec_iter /= string_vec%end() )
               string_pointer => string_vec_iter%get()
               call MAPL_GetPointer( EXPORT, ptr3d, string_pointer, _RC)
               export_clrsky = (export_clrsky .or. associated(ptr3d))
               call string_vec_iter%next()
            end do

            call string_vec%clear()
            call string_vec%push_back('OLC')
            call string_vec%push_back('OLCC5')
            call string_vec%push_back('FLNSC')
            call string_vec%push_back('LCS')
            call string_vec%push_back('LCSC5')
            string_vec_iter = string_vec%begin()
            do while ( string_vec_iter /= string_vec%end() )
               string_pointer => string_vec_iter%get()
               call MAPL_GetPointer( EXPORT, ptr2d, string_pointer, _RC)
               export_clrsky = (export_clrsky .or. associated(ptr2d))
               call string_vec_iter%next()
            end do

            ! are cloudy clean exports requested?
            export_allnoa = .false.

            call string_vec%clear()
            call string_vec%push_back('FLXA')
            call string_vec%push_back('FLXAD')
            call string_vec%push_back('FLXAU')
            string_vec_iter = string_vec%begin()
            do while ( string_vec_iter /= string_vec%end() )
               string_pointer => string_vec_iter%get()
               call MAPL_GetPointer( EXPORT, ptr3d, string_pointer, _RC)
               export_allnoa = (export_allnoa .or. associated(ptr3d))
               call string_vec_iter%next()
            end do

            call string_vec%clear()
            call string_vec%push_back('OLRA')
            call string_vec%push_back('FLNSNA')
            call string_vec%push_back('LWSA')
            string_vec_iter = string_vec%begin()
            do while ( string_vec_iter /= string_vec%end() )
               string_pointer => string_vec_iter%get()
               call MAPL_GetPointer( EXPORT, ptr2d, string_pointer, _RC)
               export_allnoa = (export_allnoa .or. associated(ptr2d))
               call string_vec_iter%next()
            end do

            ! are cloudy dirty exports requested?
            export_allsky = .false.

            call string_vec%clear()
            call string_vec%push_back('FLX')
            call string_vec%push_back('FLXD')
            call string_vec%push_back('FLXU')
            string_vec_iter = string_vec%begin()
            do while ( string_vec_iter /= string_vec%end() )
               string_pointer => string_vec_iter%get()
               call MAPL_GetPointer( EXPORT, ptr3d, string_pointer, _RC)
               export_allsky = (export_allsky .or. associated(ptr3d))
               call string_vec_iter%next()
            end do

            call string_vec%clear()
            call string_vec%push_back('OLR')
            call string_vec%push_back('SFCEM')
            call string_vec%push_back('FLNS')
            call string_vec%push_back('LWS')
            string_vec_iter = string_vec%begin()
            do while ( string_vec_iter /= string_vec%end() )
               string_pointer => string_vec_iter%get()
               call MAPL_GetPointer( EXPORT, ptr2d, string_pointer, _RC)
               export_allsky = (export_allsky .or. associated(ptr2d))
               call string_vec_iter%next()
            end do

            ! band outputs are all-sky only for the moment
            export_allsky = (export_allsky .or. any_band_output)

            ! which fluxes to calculate?
            ! the clean fluxes are also used for "dirty" fluxes if no aerosols
            calc_clrnoa = (export_clrnoa .or. (export_clrsky .and. .not.implements_aerosol_optics))
            calc_allnoa = (export_allnoa .or. (export_allsky .and. .not.implements_aerosol_optics))
            calc_clrsky =  export_clrsky
            calc_allsky =  export_allsky

            ! handle allnoa -> allsky band output when aerosols not implemented
            !   (band output currently only available for all-sky case)
            allnoa_to_allsky_band_xfer_needed = &
                 export_allsky .and. any_band_output .and. .not.implements_aerosol_optics

            ! do we actually need dirty optical properties?
            need_dirty_optical_props = &
                 (export_clrsky .or. export_allsky) .and. implements_aerosol_optics

            ! do we need cloudy optical properties?
            need_cloud_optical_props = (export_allnoa .or. export_allsky)

            ! allocation of output arrays
            if (calc_clrnoa) then
               allocate(flux_up_clrnoa(ncol,LM+1), &
                    flux_dn_clrnoa(ncol,LM+1), &
                    dfupdts_clrnoa(ncol,LM+1), _STAT)
            end if
            if (calc_allnoa) then
               allocate(flux_up_allnoa(ncol,LM+1), &
                    flux_dn_allnoa(ncol,LM+1), &
                    dfupdts_allnoa(ncol,LM+1), _STAT)
               if (allnoa_to_allsky_band_xfer_needed) then
                  allocate(bnd_flux_up_allnoa(ncol,LM+1,nbnd), &
                       bnd_dfupdts_allnoa(ncol,LM+1,nbnd), _STAT)
               end if
            end if
            if (calc_clrsky) then
               allocate(flux_up_clrsky(ncol,LM+1), &
                    flux_dn_clrsky(ncol,LM+1), &
                    dfupdts_clrsky(ncol,LM+1), _STAT)
            end if
            if (calc_allsky) then
               allocate(flux_up_allsky(ncol,LM+1), &
                    flux_dn_allsky(ncol,LM+1), &
                    dfupdts_allsky(ncol,LM+1), _STAT)
               if (any_band_output) then
                  allocate(bnd_flux_up_allsky(ncol,LM+1,nbnd), &
                       bnd_dfupdts_allsky(ncol,LM+1,nbnd), _STAT)
               end if
            end if

            ! IMPORTANT: Specify the type (#streams) of the LW RT calculations in clean_optical_props
            ! While the cloud optics file currently provides two-stream properties, as does the
            ! aerosol system, we may choose any number of streams for the actual RT calculations by
            ! the appropriate instantiation of clean_optical_props here. The increment() statements
            ! below implicitly convert all component optical properties to this number of streams.
            ! Everything else in the code should adapt polymorphically without modification.
            ! options are: 1scl (no scattering), 2str (2-stream), or nstr (n-stream)
            ! For 1scl, must also specify the number of Gauss angles (nga) below.
            ! For nstr, must also specify the number of phase function moments (nmom) below.
            ! After Feb2020 update:
            ! Even for optical_props_2str, the default rte method is to use rescaled LW transport
            !   to account for scattering (in which case nga is used). To force explicit 2-stream
            !   scattering, must select u2s = .true. and allocate optical_props_2str below.

            ! LW uses ty_optical_props_2str (see PROCESS_RRTMGP_LW_BLOCK).
            ! nga, nmom, u2s are determined here once and passed to each block call.
            nga  = 1    ! used when not u2s; must be >= 1
            nmom = 2    ! used only if nstr; must be >= 2
            u2s  = .false.
            call MAPL_GetResource( &
                 MAPL, u2s ,'RRTMGP_LW_USE_2STREAM:',    DEFAULT=u2s, _RC)
            _ASSERT(.not.u2s,'lw_solver_2stream() does not currently support Jacobians')
            call MAPL_GetResource( &
                 MAPL, nga ,'RRTMGP_LW_N_GAUSS_ANGLES:', DEFAULT=nga, _RC)

            ! get cloud optical properties (band-only)
            if (need_cloud_optical_props) then
               ! pmn: some of this should be done only once per run

               ! load and init cloud_optics from file:
               ! gets appropriate coefficients needed to calculate
               ! cloud optical properties from cloud physical properties
               call MAPL_GetResource( &
                    MAPL, cloud_optics_file, "RRTMGP_CLOUD_OPTICS_LW:", &
                    DEFAULT='rrtmgp-clouds-lw.nc', _RC)
               call MAPL_GetResource( &
                    MAPL, cloud_optics_type, "RRTMGP_CLOUD_OPTICS_TYPE_LW:", &
                    DEFAULT='LUT', _RC)
               call MAPL_TimerOn(MAPL,"---RRTMGP_IO_CLOUDS",_RC)
               if (trim(cloud_optics_type)=='LUT') then
                  call load_cld_lutcoeff (cloud_optics, cloud_optics_file)
               elseif (trim(cloud_optics_type)=='PADE') then
                  call load_cld_padecoeff(cloud_optics, cloud_optics_file)
               else
                  TEST_('unknown cloud_optics_type: '//trim(cloud_optics_file))
               end if
               call MAPL_TimerOff(MAPL,"---RRTMGP_IO_CLOUDS",_RC)

               ! ice surface roughness category for Yang (2013) ice optics
               ! icergh: 1 = none, 2 = medium, 3 = high
               call MAPL_GetResource( &
                    MAPL, icergh, "RRTMGP_ICE_ROUGHNESS_LW:", &
                    DEFAULT=2, _RC)
               TEST_(cloud_optics%set_ice_roughness(icergh))

               ! cloud optics file is currently two-stream; cloud_props_bnd/gpt
               ! are allocated/init'd per-block inside PROCESS_RRTMGP_LW_BLOCK.

               ! read desired cloud overlap type

               call MAPL_GetResource( &
                    MAPL, cloud_overlap_type, "RRTMGP_CLOUD_OVERLAP_TYPE_LW:", &
                    DEFAULT='GEN_MAX_RAN_OVERLAP', _RC)

               ! GEN_MAX_RAN_OVERLAP uses correlation lengths
               !   and possibly inhomogeneous condensate
               gen_mro = (cloud_overlap_type == "GEN_MAX_RAN_OVERLAP")
               if (gen_mro) then

                  ! condensate inhomogeneous?
                  ! see RadiationGC initialization
                  cond_inhomo = condensate_inhomogeneous()

                  ! Compute decorrelation length scales [m]
                  allocate(adl(ncol),_STAT)
                  call correlation_length_cloud_fraction(ncol, ncol, doy, reshape(LATS,(/ncol/)), adl)
                  if (cond_inhomo) then
                     allocate(rdl(ncol),_STAT)
                     call correlation_length_condensate(ncol, ncol, doy, reshape(LATS,(/ncol/)), rdl)
                  endif

               endif

               ! Random number setup:
               !   We will use the Philox4x32-10 or ARS5 BRNGs from MKL VSL.
               !   Both are keyed families of counter-based PRNGs with a large period 2^130
               ! and a minimal state space (unlike the large Mersenne Twister state).
               !   Philox4x32-10 has a 64-bit key and a 128-bit counter and is very fast on GPUs.
               !   ARS5 has a 128-bit key and a 128-bit counter and is superfast on CPUs for
               ! which AES-NI instructions are hardware implemented.
               !   The SEEDING STRATEGY we will follow is to use a unique key for the gricolumn
               ! location and the simulation time. This gives a repeatable set of random numbers
               ! that remains the same for the members of an ensemble. If a different set is
               ! required for ensemble members, then the model state, such as the fractional
               ! part of the surface pressure, should be incorporated into the key.
               !   To get a different set of random numbers for the SW, for example, either a
               ! key change or a counter advance will be needed.
               !
               ! Time Component of key:
               ! ~~~~~~~~~~~~~~~~~~~~~~
               ! 1. No need to update more frequently than once per LW refresh.
               ! 2. should reference the number of such intervals since a fixed time, so
               !   that agnostic to stop/restart schedule.
               !
               ! Space component of key:
               ! ~~~~~~~~~~~~~~~~~~~~~~~
               ! 1. should be based on some globally unique index for a gridcolumn, so that
               !   each gridcolumn is independent and so it is agnostic to runs with varying
               !   decompositions among processors.
               ! 2. 2^32 = 4,294,967,296 or about 2.1475e9 positives, which can represent
               !   globe at over 1/180th degree resolution, so plenty for forseeable
               !   future.
               !
               ! Philox seeding:
               ! ~~~~~~~~~~~~~~~
               !   1. a scalar 32-bit seed sets the lower bits of the key k.
               !   2. a vector of 32-bit seeds of length N is used as follows:
               !      (a) N = 0:      k = c = 0;
               !      (b) N in {1,2}: seeds(1(:2)) set lower (and upper) words of key
               !      (c) N > 2:      ditto plus seeds(3:min(N,6)) set counter c,
               !                        starting from lowest word and working up.
               !
               ! Estimate of maximum LW random numbers needed per gridcolumn:
               ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               ! For the current clouds, with homogeneous optical properies in the cloudy part
               ! of each layer, only one random number per layer and gpt is needed. So estimate
               ! LM * ngpt <~ 132 * 256 = 33,792 < 2^16 = 65,536

               allocate(seeds(3),_STAT) ! 2-word key plus word1 of counter

               ! seed(1), the column part (word1) of key is set later
               ! but get required global indicies of local rectangular grid here
               call MAPL_GridGet(ESMFGRID, globalCellCountPerDim=Gdims, _RC)
               IM_World = Gdims(1); JM_World = Gdims(2)
               call MAPL_GridGetInterior (ESMFGRID,iBeg,iEnd,jBeg,jEnd)

               ! get time part (word2) of key
               call ESMF_ClockGet(CLOCK, currTIME=CurrentTime, _RC)
               call ESMF_TimeSet (ReferenceTime, yy=2000, mm=1, dd=1, _RC)
               call ESMF_AlarmGet(ALARM, RINGINTERVAL=RefreshInterval, _RC)
               seeds(2) = int((CurrentTime - ReferenceTime) / RefreshInterval)

               ! for LW start at counter=0
               seeds(3) = 0

               ! get a view of cloud inputs with collapsed horizontal dimensions
               call c_f_pointer(c_loc(CWC), CWC_3d, [IM*JM,LM,5])
               call c_f_pointer(c_loc(REFF),REFF_3d,[IM*JM,LM,5])

            end if ! need_cloud_optical_props

            ! set aerosol optical properties
            if (need_dirty_optical_props) then

               ! aerosol optics system is currently two-stream
               ! aer_props alloc+init handled per-block in PROCESS_RRTMGP_LW_BLOCK.
               ! get a view of aerosol system inputs with collapsed horizontal dimensions
               ! (aer_props is always ty_optical_props_2str for LW)
               call c_f_pointer(c_loc(TAUA),TAUA_3d,[IM*JM,LM,NB_IRRAD])
               call c_f_pointer(c_loc(SSAA),SSAA_3d,[IM*JM,LM,NB_IRRAD])
               call c_f_pointer(c_loc(ASYA),ASYA_3d,[IM*JM,LM,NB_IRRAD])

            end if

            !-------------------------------------------------------!
            ! Loop over blocks of blockSize columns                 !
            !  - choose rrtmgp_blockSize for memory/time efficiency !
            !  - all blocks including the final partial block are   !
            !    handled uniformly using ceiling division           !
            !-------------------------------------------------------!

            call MAPL_GetResource( MAPL, &
                 rrtmgp_blockSize, "RRTMGP_LW_BLOCKSIZE:", DEFAULT=4, _RC)
            _ASSERT(rrtmgp_blockSize >= 1,'invalid RRTMGP_LW_BLOCKSIZE')

            ! Total number of blocks including any final partial block
            nBlocks = (ncol + rrtmgp_blockSize - 1) / rrtmgp_blockSize

            ! loop over all blocks
            loop_status = ESMF_SUCCESS
            !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(STATUS)
            do b = 1, nBlocks
               call PROCESS_RRTMGP_LW_BLOCK( &
                    b, rrtmgp_blockSize, ncol, LM, nmom, ngpt, nga, &
                    IM, IM_World, iBeg, jBeg, &
                    top_at_1, u2s, &
                    seeds(2), seeds(3), &
                    cwp_fac, &
                    need_cloud_optical_props, need_dirty_optical_props, &
                    gen_mro, cond_inhomo, cloud_overlap_type, &
                    calc_clrnoa, calc_allnoa, calc_clrsky, calc_allsky, &
                    allnoa_to_allsky_band_xfer_needed, any_band_output, &
                    export_clrsky, export_allsky, implements_aerosol_optics, &
                    k_dist, cloud_optics, gas_concs, &
                    p_lay, p_lev, t_lay, t_lev, t_sfc, dp_wp, cf_wp, dzmid, emis_sfc, &
                    adl=adl, rdl=rdl, &
                    CWC_3d=CWC_3d, REFF_3d=REFF_3d, &
                    TAUA_3d=TAUA_3d, SSAA_3d=SSAA_3d, ASYA_3d=ASYA_3d, &
                    flux_up_clrnoa=flux_up_clrnoa, &
                    flux_dn_clrnoa=flux_dn_clrnoa, &
                    dfupdts_clrnoa=dfupdts_clrnoa, &
                    flux_up_allnoa=flux_up_allnoa, &
                    flux_dn_allnoa=flux_dn_allnoa, &
                    dfupdts_allnoa=dfupdts_allnoa, &
                    bnd_flux_up_allnoa=bnd_flux_up_allnoa, &
                    bnd_dfupdts_allnoa=bnd_dfupdts_allnoa, &
                    flux_up_clrsky=flux_up_clrsky, &
                    flux_dn_clrsky=flux_dn_clrsky, &
                    dfupdts_clrsky=dfupdts_clrsky, &
                    flux_up_allsky=flux_up_allsky, &
                    flux_dn_allsky=flux_dn_allsky, &
                    dfupdts_allsky=dfupdts_allsky, &
                    bnd_flux_up_allsky=bnd_flux_up_allsky, &
                    bnd_dfupdts_allsky=bnd_dfupdts_allsky, &
                    MAPL=MAPL, RC=STATUS)
               if (STATUS /= ESMF_SUCCESS) then
                  !$OMP ATOMIC WRITE
                  loop_status = STATUS
               end if
            end do ! loop over blocks
            !$OMP END PARALLEL DO
            VERIFY_(loop_status)

            ! tidy up
            if (need_dirty_optical_props) nullify(TAUA_3d,SSAA_3d,ASYA_3d)
            if (need_cloud_optical_props) nullify(CWC_3d,REFF_3d)

            call MAPL_TimerOn(MAPL,"---RRTMGP_POST",_RC)

            ! load output arrays
            ! note: the upward fluxes must be NEGATED for the downward +ve conventionS
            ! likewise, the DFDTS* are the derivatives of the NEGATED upward fluxes wrt TS
            if (export_clrnoa) then
               FLAU_INT = real(reshape(-flux_up_clrnoa, (/IM,JM,LM+1/)))
               FLAD_INT = real(reshape( flux_dn_clrnoa, (/IM,JM,LM+1/)))
               DFDTSCNA = real(reshape(-dfupdts_clrnoa, (/IM,JM,LM+1/)))
            end if
            if (export_allnoa) then
               FLXAU_INT = real(reshape(-flux_up_allnoa, (/IM,JM,LM+1/)))
               FLXAD_INT = real(reshape( flux_dn_allnoa, (/IM,JM,LM+1/)))
               DFDTSNA   = real(reshape(-dfupdts_allnoa, (/IM,JM,LM+1/)))
            end if
            if (export_clrsky) then
               FLCU_INT = real(reshape(-flux_up_clrsky, (/IM,JM,LM+1/)))
               FLCD_INT = real(reshape( flux_dn_clrsky, (/IM,JM,LM+1/)))
               DFDTSC   = real(reshape(-dfupdts_clrsky, (/IM,JM,LM+1/)))
            end if
            if (export_allsky) then
               FLXU_INT = real(reshape(-flux_up_allsky, (/IM,JM,LM+1/)))
               FLXD_INT = real(reshape( flux_dn_allsky, (/IM,JM,LM+1/)))
               DFDTS    = real(reshape(-dfupdts_allsky, (/IM,JM,LM+1/)))
            end if

            !mjs: Corrected emitted at the surface to remove reflected
            !     from upward. Note that emiss is the same for all bands,
            !     so we use band 1 for the total flux.
            SFCEM_INT = real(reshape( &
                 -flux_up_allsky(:,LM+1) + flux_dn_allsky(:,LM+1) * (1._wp - emis_sfc(1,:)), &
                 (/IM,JM/)))

            ! band OLR and Tsfc Jacobian
            ! These are direct INTERNALs and do not need the above negation for upward fluxes
            if (export_allsky) then
               do ib = 1,nbnd
                  if (band_output(ib)) then
                     write(bb,'(I0.2)') ib
                     call MAPL_GetPointer(INTERNAL, ptr2d, 'OLRB'//bb//'RG', _RC)
                     ptr2d = real(reshape(bnd_flux_up_allsky(:,1,ib), [IM,JM]))
                     call MAPL_GetPointer(INTERNAL, ptr2d, 'DOLRB'//bb//'RGDT', _RC)
                     ptr2d = real(reshape(bnd_dfupdts_allsky(:,1,ib), [IM,JM]))
                  end if
               end do
            end if

            ! clean up
            deallocate(t_sfc,emis_sfc,_STAT)
            deallocate(p_lay,t_lay,p_lev,t_lev,dp_wp,cf_wp,dzmid,_STAT)
            if (need_cloud_optical_props) then
               deallocate(seeds,_STAT)
               if (gen_mro) then
                  deallocate(adl,_STAT)
                  if (cond_inhomo) then
                     deallocate(rdl,_STAT)
                  endif
               endif
               call cloud_optics%finalize()
            end if
            if (calc_clrnoa) then
               deallocate(flux_up_clrnoa, flux_dn_clrnoa, dfupdts_clrnoa, _STAT)
            end if
            if (calc_allnoa) then
               deallocate(flux_up_allnoa, flux_dn_allnoa, dfupdts_allnoa, _STAT)
               if (allnoa_to_allsky_band_xfer_needed) then
                  deallocate(bnd_flux_up_allnoa, bnd_dfupdts_allnoa, _STAT)
               end if
            end if
            if (calc_clrsky) then
               deallocate(flux_up_clrsky, flux_dn_clrsky, dfupdts_clrsky, _STAT)
            end if
            if (calc_allsky) then
               deallocate(flux_up_allsky, flux_dn_allsky, dfupdts_allsky, _STAT)
               if (any_band_output) then
                  deallocate(bnd_flux_up_allsky, bnd_dfupdts_allsky, _STAT)
               end if
            end if

            call MAPL_TimerOff(MAPL,"---RRTMGP_POST",_RC)

            call MAPL_TimerOff(MAPL,"--RRTMGP",_RC)

         else if (USE_RRTMG) then

            call MAPL_TimerOn(MAPL,"--RRTMG",_RC)

            if (LM > 72) then
               call MAPL_GetResource(MAPL,USE_PRECIP_IN_RADIATION,'RRTMGLW_USE_PRECIP_IN_RADIATION:',DEFAULT=.TRUE.,_RC)
            else
               call MAPL_GetResource(MAPL,USE_PRECIP_IN_RADIATION,'RRTMGLW_USE_PRECIP_IN_RADIATION:',DEFAULT=.FALSE.,_RC)
            endif

            call MAPL_GetResource(MAPL,PARTITION_SIZE,'RRTMGLW_PARTITION_SIZE:',DEFAULT=4,_RC)

            ! reversed profiles for RRTMG (1=bottom layer)
            ! note 0:LM indexing for [PT]LEV_R
            ! but 1:LM+1 for [UD]FLX[C] and DUFLX[C]_DTS
            allocate(FCLD_R(IM*JM,LM),_STAT)
            allocate(TLEV_R(IM*JM,0:LM),_STAT)
            allocate(PLE_R(IM*JM,0:LM),_STAT)
            allocate(ZM_R(IM*JM,LM),_STAT)
            allocate(EMISS(IM*JM,NB_RRTMG),_STAT)
            allocate(CLIQWP(IM*JM,LM),_STAT)
            allocate(CICEWP(IM*JM,LM),_STAT)
            allocate(RELIQ(IM*JM,LM),_STAT)
            allocate(REICE(IM*JM,LM),_STAT)
            allocate(TAUAER(IM*JM,LM,NB_RRTMG),_STAT)
            allocate(PL_R(IM*JM,LM),_STAT)
            allocate(T_R(IM*JM,LM),_STAT)
            allocate(Q_R(IM*JM,LM),_STAT)
            allocate(O2_R(IM*JM,LM),_STAT)
            allocate(O3_R(IM*JM,LM),_STAT)
            allocate(CO2_R(IM*JM,LM),_STAT)
            allocate(CH4_R(IM*JM,LM),_STAT)
            allocate(N2O_R(IM*JM,LM),_STAT)
            allocate(CFC11_R(IM*JM,LM),_STAT)
            allocate(CFC12_R(IM*JM,LM),_STAT)
            allocate(CFC22_R(IM*JM,LM),_STAT)
            allocate(CCL4_R(IM*JM,LM),_STAT)
            allocate(TSFC(IM*JM),_STAT)
            allocate(UFLX(IM*JM,LM+1),_STAT)
            allocate(DFLX(IM*JM,LM+1),_STAT)
            allocate(UFLXC(IM*JM,LM+1),_STAT)
            allocate(DFLXC(IM*JM,LM+1),_STAT)
            allocate(DUFLX_DTS(IM*JM,LM+1),_STAT)
            allocate(DUFLXC_DTS(IM*JM,LM+1),_STAT)
            allocate(CLEARCOUNTS(IM*JM,4),_STAT)
            allocate(ALAT(IM*JM),_STAT)
            allocate(OLRBRG      (nbndlw,IM*JM),_STAT)
            allocate(DOLRBRG_DTS (nbndlw,IM*JM),_STAT)

            ! choices for cloud physical to optical conversion
            call MAPL_GetResource(MAPL,ICEFLGLW,'RRTMG_ICEFLG:',DEFAULT=3,_RC)
            call MAPL_GetResource(MAPL,LIQFLGLW,'RRTMG_LIQFLG:',DEFAULT=1,_RC)

            ! calculate derivatives of upward flux with Tsurf
            Ts_derivs = .true.

            call MAPL_TimerOn(MAPL,"---RRTMG_FLIP",_RC)

            ! reverse super-layer interface indicies
            LCLDMH = LM - LCLDMH + 1
            LCLDLM = LM - LCLDLM + 1

            ! collapse horizontal indicies and flip in vertical
            !   (RRTMG indexed bottom to top)
            IJ = 0
            do J = 1,JM
            do I = 1,IM
               IJ = IJ + 1

               TSFC (IJ)   = TS  (I,J)
               EMISS(IJ,:) = EMIS(I,J) ! all bands get same emissivity
               ALAT (IJ)   = LATS(I,J)

               ! calculation of level temperature (still in model ordering)
               ! note: PLE(0:LM) but TLEV(1:LM+1)
               DP(1) = PLE(I,J,1)-PLE(I,J,0)
               do K = 2,LM
                  DP(K) = (PLE(I,J,K)-PLE(I,J,K-1) )
                  TLEV(K) = (T(I,J,K-1) * DP(K) + T(I,J,K) * DP(K-1)) &
                       / (DP(K-1) + DP(K))
               enddo
               TLEV(LM+1) =  TS(I,J) ! 'surface'
               TLEV(   1) = TLEV(2)  ! model top

               !  Flip in vertical
               do K = 1,LM
                  LV = LM-K+1  ! LM --> 1

                  ! Convert content [kg/kg] to path [g/m2]
                  ! using hydrostatic eqn dp/g ~ rho*dz,
                  ! so conversion factor is 1000*dp/g ~ 1.02*100*dp.
                  ! pmn: why not use MAPL_GRAV explicitly?
                  xx = 1.02*100*DP(LV)
                  if (USE_PRECIP_IN_RADIATION) then
                     LWT = CWC(I,J,LV,KLIQUID)+CWC(I,J,LV,KRAIN)
                     CLIQWP(IJ,K) = xx*(LWT)
                     if (LWT > 0.0) then
                        RELIQ (IJ,K) = ( REFF(I,J,LV,KLIQUID)*CWC(I,J,LV,KLIQUID) + &
                             REFF(I,J,LV,KRAIN  )*CWC(I,J,LV,KRAIN  ) ) / LWT
                     else
                        RELIQ (IJ,K) = 14.0
                     endif
                     IWT = CWC(I,J,LV,KICE)+CWC(I,J,LV,KSNOW)+CWC(I,J,LV,KGRAUPEL)
                     CICEWP(IJ,K) = xx*(IWT)
                     if (IWT > 0.0) then
                        REICE (IJ,K) = ( REFF(I,J,LV,KICE    )*CWC(I,J,LV,KICE    ) + &
                             REFF(I,J,LV,KSNOW   )*CWC(I,J,LV,KSNOW   ) + &
                             REFF(I,J,LV,KGRAUPEL)*CWC(I,J,LV,KGRAUPEL) ) / IWT
                     else
                        REICE (IJ,K) = 36.0
                     endif
                  else
                     CLIQWP(IJ,K) = xx*CWC(I,J,LV,KLIQUID)
                     CICEWP(IJ,K) = xx*CWC(I,J,LV,KICE)
                     RELIQ (IJ,K) =   REFF(I,J,LV,KLIQUID)
                     REICE (IJ,K) =   REFF(I,J,LV,KICE   )
                  endif

                  ! impose RRTMG re_liq limits
                  if    (LIQFLGLW.eq.0) then
                     ! pmn: this one not available inside RRTMG_LW
                     RELIQ(IJ,K) = min(max(RELIQ(IJ,K),5.0),10.0)
                  elseif (LIQFLGLW.eq.1) then
                     RELIQ(IJ,K) = min(max(RELIQ(IJ,K),2.5),60.0)
                  endif

                  ! impose RRTMG re_ice limits
                  if     (ICEFLGLW.eq.0) then
                     REICE(IJ,K) = min(max(REICE(IJ,K),10.0),30.0)
                  elseif (ICEFLGLW.eq.1) then
                     REICE(IJ,K) = min(max(REICE(IJ,K),13.0),130.0)
                  elseif (ICEFLGLW.eq.2) then
                     REICE(IJ,K) = min(max(REICE(IJ,K), 5.0),131.0)
                  elseif (ICEFLGLW.eq.3) then
                     REICE(IJ,K) = min(max(REICE(IJ,K), 5.0),140.0)
                  elseif (ICEFLGLW.eq.4) then
                     REICE(IJ,K) = min(max(REICE(IJ,K)*2.,1.0),200.0)
                  endif

                  ! flipping for LEVEL quantities
                  ! PLE_R(0:LM) = PLE(LM:0)
                  ! TLEV_R(0:LM) = TLEV(LM+1:1)
                  ! top-of-model LEVEL (RRTMG LM) done later
                  PLE_R  (IJ,K-1) = PLE(I,J,LV)/100. ! [hPa]
                  TLEV_R (IJ,K-1) = TLEV(LV+1)

                  ! more flipping for layer quantities
                  ! Q [specific humidity] --> Q_R [volume mixing ratio]
                  ! O3 [mass mixing ratio] --> O3_R [volume mixing ratio]
                  PL_R   (IJ,K) = PL(I,J,LV)/100.  ! [hPa]
                  T_R    (IJ,K) = T(I,J,LV)
                  Q_R    (IJ,K) = Q(I,J,LV) / (1.-Q(I,J,LV)) * (MAPL_AIRMW/MAPL_H2OMW)
                  O3_R   (IJ,K) = O3(I,J,LV) * (MAPL_AIRMW/MAPL_O3MW)
                  CH4_R  (IJ,K) = CH4(I,J,LV)
                  N2O_R  (IJ,K) = N2O(I,J,LV)
                  if (associated(CO2_3d)) then ! <<>> MSL
                     CO2_R  (IJ,k) = CO2_3d(I,J,LV)
                  else
                     CO2_R  (IJ,k) = CO2_FIXED
                  endif
                  O2_R   (IJ,K) = O2
                  CCL4_R (IJ,K) = CCL4
                  CFC11_R(IJ,K) = CFC11(I,J,LV)
                  CFC12_R(IJ,K) = CFC12(I,J,LV)
                  CFC22_R(IJ,K) = HCFC22(I,J,LV)
                  FCLD_R (IJ,K) = FCLD(I,J,LV)

                  ! RRTMG_LW does not scatter, so pass ABSORPTION aerosol
                  ! optical thickness to RRTMG. Remember that SSAA is the
                  ! aerosol system's *un*-normalized single scattering albedo,
                  ! which is actually tau_ext * omega0 = tau_scat, and TAUA
                  ! is the aerosol extinction optical thickness.
                  ! PMN 2022-01-19 Added max(,0.) ... it shouldn't happen
                  !   that tau_ext < tau_scat, but since these TAUA and SSAA
                  !   come directly from the radiatively active aerosols
                  !   system, we provide a simple mitigation here.
                  TAUAER(IJ,K,:) = max(TAUA(I,J,LV,:) - SSAA(I,J,LV,:), 0.)

               enddo

               ! finish off top-of-model LEVEL
               PLE_R (IJ,LM) = PLE(I,J,0)/100. ! [hPa]
               TLEV_R(IJ,LM) = TLEV(1)

               ! Calculate the LAYER (mid-point) heights.
               ! The interlayer distances are needed for the calculations
               ! of inter-layer correlation for cloud overlapping in RRTMG.
               ! Only *relative* distances matter, so wolog set ZM_R(1) = 0.
               ! pmn: 2021-04-21 this calculation was wrong in earlier revisions.
               ZM_R(IJ,1) = 0.
               do K=2,LM
                  ! dz ~ RT/g x dp/p by hysrostatic eqn and ideal gas eqn.
                  ! The jump from LAYER k-1 to k is centered on LEVEL k-1
                  !   since the RRTMG LEVEL (LE[V]_R) indices are zero-based
                  ZM_R(IJ,K) = ZM_R(IJ,K-1) + MAPL_RGAS*TLEV_R(IJ,K-1)/MAPL_GRAV &
                       * (PL_R(IJ,K-1)-PL_R(IJ,K))/PLE_R(IJ,K-1)
               enddo

            enddo ! IM
            enddo ! JM

            ! Clean up negatives
            WHERE (Q_R < 0.) Q_R = 0.
            WHERE (O3_R < 0.) O3_R = 0.
            WHERE (CH4_R < 0.) CH4_R = 0.
            WHERE (N2O_R < 0.) N2O_R = 0.
            WHERE (CO2_R < 0.) CO2_R = 0.
            WHERE (O2_R < 0.) O2_R = 0.
            WHERE (CCL4_R < 0.) CCL4_R = 0.
            WHERE (CFC11_R < 0.) CFC11_R = 0.
            WHERE (CFC12_R < 0.) CFC12_R = 0.
            WHERE (CFC22_R < 0.) CFC22_R = 0.
            WHERE (FCLD_R < 0.) FCLD_R = 0.

            call MAPL_TimerOff(MAPL,"---RRTMG_FLIP",_RC)

            call MAPL_TimerOn(MAPL,"---RRTMG_INIT",_RC)

            ! pmn: consider putting futher up calling tree?
            ! pmn: only needs to be done once per run, but does consume memory
            call RRTMG_LW_INI

            call MAPL_TimerOff(MAPL,"---RRTMG_INIT",_RC)

            call MAPL_TimerOn(MAPL,"---RRTMG_RUN",_RC)

            if (nRATS .gt. 0) then !<<>> MSL
               allocate(UFLXRAT(IM*JM,LM+1,nRATS),      _STAT)
               allocate(DFLXRAT(IM*JM,LM+1,nRATS),      _STAT)
               allocate(DUFLX_DT_RAT(IM*JM,LM+1,nRATS), _STAT)
               !         allocate(DFDTS_RAT(IM,JM,LM+1,nRATS),    _STAT)
               !         allocate(FLXU_INT_RAT(IM,JM,LM+1,nRATS),    _STAT)
               !         allocate(FLXD_INT_RAT(IM,JM,LM+1,nRATS),    _STAT)
               !         allocate(FLX_INT_RAT(IM,JM,LM+1,nRATS) ,    _STAT)
               call MAPL_GetPointer(INTERNAL, DFDTS_RAT,     'DFDTS_RAT',  RC=STATUS)
               call MAPL_GetPointer(INTERNAL, FLX_INT_RAT,   'FLX_RAT',    RC=STATUS)
               call MAPL_GetPointer(INTERNAL, FLXU_INT_RAT,  'FLXU_RAT',   RC=STATUS)
               call MAPL_GetPointer(INTERNAL, FLXD_INT_RAT,  'FLXD_RAT',   RC=STATUS)
               call MAPL_GetPointer(INTERNAL, SFCEM_INT_RAT, 'SFCEM_RAT',  RC=STATUS)
            endif

            ! Begin analysis for RATS-specific rad diagnostics
            do n = 1,nRATS !<<>> MSL

               ! Zero out the correct RAT gas field
               ! O3 CO2 CH4 N2O CFC11 CFC12 CFC22 CCl4
               select case(nameRATS(n))
               case('H2O')
                  TMP_R = Q_R
                  Q_R = 0.e0
               case('O3')
                  TMP_R = O3_R
                  O3_R = 0.e0
               case('CO2')
                  TMP_R = CO2_R
                  !            CO2_R = CO2_FIXED! Testing
                  CO2_R = 0.e0
                  !            CO2_R = CO2_R*0.99e0
               case('CH4')
                  TMP_R = CH4_R
                  CH4_R = 0.e0
               case('N2O')
                  TMP_R = N2O_R
                  N2O_R = 0.e0
               case('CFC11')
                  TMP_R = CFC11_R
                  CFC11_R = 0.e0
               case('CFC12')
                  TMP_R = CFC12_R
                  CFC12_R = 0.e0
               case('HCFC22_R')
                  TMP_R = CFC22_R
                  CFC22_R = 0.e0
               end select

               ! Call the long wave code with a given RAT set to zero
               call RRTMG_LW (IM*JM, LM, PARTITION_SIZE, TS_DERIVS, &
                    PL_R, PLE_R, T_R, TLEV_R, TSFC, EMISS, &
                    Q_R, O3_R, CO2_R, CH4_R, N2O_R, O2_R, &
                    CFC11_R, CFC12_R, CFC22_R, CCL4_R, &
                    FCLD_R, CICEWP, CLIQWP, REICE, RELIQ, ICEFLGLW, LIQFLGLW, &
                    TAUAER, ZM_R, ALAT, DOY, LCLDLM, LCLDMH, CLEARCOUNTS, &
                    UFLXRAT(:,:,n), DFLXRAT(:,:,n), UFLXC, DFLXC, DUFLX_DT_RAT(:,:,n), DUFLXC_DTS, &
                    BAND_OUTPUT, OLRBRG, DOLRBRG_DTS)

               ! Make sure to set the RAT gas column back to correct vals
               select case(nameRATS(n))
               case('H2O')
                  Q_R = TMP_R
               case('O3')
                  O3_R = TMP_R
               case('CO2')
                  CO2_R = TMP_R
               case('CH4')
                  CH4_R = TMP_R
               case('N2O')
                  N2O_R = TMP_R
               case('CFC11')
                  CFC11_R = TMP_R
               case('CFC12')
                  CFC12_R = TMP_R
               case('HCFC22')
                  CFC22_R = TMP_R
               end select

            enddo
            ! <<>> end RATS analysis

            call RRTMG_LW (IM*JM, LM, PARTITION_SIZE, TS_DERIVS, &
                 PL_R, PLE_R, T_R, TLEV_R, TSFC, EMISS, &
                 Q_R, O3_R, CO2_R, CH4_R, N2O_R, O2_R, &
                 CFC11_R, CFC12_R, CFC22_R, CCL4_R, &
                 FCLD_R, CICEWP, CLIQWP, REICE, RELIQ, ICEFLGLW, LIQFLGLW, &
                 TAUAER, ZM_R, ALAT, DOY, LCLDLM, LCLDMH, CLEARCOUNTS, &
                 UFLX, DFLX, UFLXC, DFLXC, DUFLX_DTS, DUFLXC_DTS, &
                 BAND_OUTPUT, OLRBRG, DOLRBRG_DTS)

            call MAPL_TimerOff(MAPL,"---RRTMG_RUN",_RC)

            call MAPL_TimerOn(MAPL,"---RRTMG_FLIP",_RC)

            ! for outputs, unpack flattened horizontal and flip back vertical
            IJ = 0
            do J = 1,JM
            do I = 1,IM
               IJ = IJ + 1

               ! convert super-layer clearCounts to cloud fractions
               if(associated(CLDTTLW)) then
                  CLDTTLW(I,J) = 1.0 - CLEARCOUNTS(IJ,1)/float(NGPTLW)
               endif
               if(associated(CLDHILW)) then
                  CLDHILW(I,J) = 1.0 - CLEARCOUNTS(IJ,2)/float(NGPTLW)
               endif
               if(associated(CLDMDLW)) then
                  CLDMDLW(I,J) = 1.0 - CLEARCOUNTS(IJ,3)/float(NGPTLW)
               endif
               if(associated(CLDLOLW)) then
                  CLDLOLW(I,J) = 1.0 - CLEARCOUNTS(IJ,4)/float(NGPTLW)
               endif

               ! upward negative in GEOS-5 convention
               do K = 0,LM
                  LV = LM-K+1
                  FLXU_INT(I,J,K) =-UFLX      (IJ,LV)
                  FLXD_INT(I,J,K) = DFLX      (IJ,LV)
                  FLCU_INT(I,J,K) =-UFLXC     (IJ,LV)
                  FLCD_INT(I,J,K) = DFLXC     (IJ,LV)
                  DFDTS   (I,J,K) =-DUFLX_DTS (IJ,LV)
                  DFDTSC  (I,J,K) =-DUFLXC_DTS(IJ,LV)
               enddo

               ! Reflected LW is not counted in surface emitted. Also, for now,
               ! surface emitted is positive downwards consistent with Chou-Suarez.
               ! (Note: All bands use the same emissivity)
               SFCEM_INT(I,J) = -( UFLX(IJ,1) - DFLX(IJ,1)*(1.-EMIS(I,J)) )

               if (nRATS .gt. 0) then !<<>> MSL
                  do k=0,LM
                     LV = LM-k+1
                     FLXU_INT_RAT(i,j,k,:) =-UFLXRAT     (IJ,LV,:)
                     FLXD_INT_RAT(i,j,k,:) = DFLXRAT     (IJ,LV,:)
                     DFDTS_RAT   (i,j,k,:) =-DUFLX_DT_RAT(IJ,LV,:)
                  enddo
                  SFCEM_INT_RAT(i,j,:) = UFLXRAT(IJ,1,:) - DFLXRAT(IJ,1,:)*(1.0-EMISS(IJ,1))
               endif

            enddo ! IM
            enddo ! JM

            ! band OLR and brightness temperatures
            do ibnd = 1,nbndlw
               if (band_output(ibnd)) then
                  write(bb,'(I0.2)') ibnd

                  call MAPL_GetPointer(INTERNAL, ptr2d, 'OLRB'//bb//'RG', _RC)
                  ptr2d = reshape(OLRBRG (ibnd,:), [IM,JM])

                  call MAPL_GetPointer(INTERNAL, ptr2d, 'DOLRB'//bb//'RGDT', _RC)
                  ptr2d = reshape(DOLRBRG_DTS (ibnd,:), [IM,JM])

               end if
            end do

            call MAPL_TimerOff(MAPL,"---RRTMG_FLIP",_RC)

            ! pmn:
            ! RRTMG does not provide no-aerosol derivatives
            ! so set no-aerosol to aerosol derivatives
            DFDTSNA  = DFDTS
            DFDTSCNA = DFDTSC

            deallocate(FCLD_R,_STAT)
            deallocate(TLEV_R,_STAT)
            deallocate(PLE_R,_STAT)
            deallocate(ZM_R,_STAT)
            deallocate(EMISS,_STAT)
            deallocate(CLIQWP,_STAT)
            deallocate(CICEWP,_STAT)
            deallocate(RELIQ,_STAT)
            deallocate(REICE,_STAT)
            deallocate(TAUAER,_STAT)
            deallocate(PL_R,_STAT)
            deallocate(T_R,_STAT)
            deallocate(Q_R,_STAT)
            deallocate(O2_R,_STAT)
            deallocate(O3_R,_STAT)
            deallocate(CO2_R,_STAT)
            deallocate(CH4_R,_STAT)
            deallocate(N2O_R,_STAT)
            deallocate(CFC11_R,_STAT)
            deallocate(CFC12_R,_STAT)
            deallocate(CFC22_R,_STAT)
            deallocate(CCL4_R,_STAT)
            deallocate(TSFC,_STAT)
            deallocate(UFLX,_STAT)
            deallocate(DFLX,_STAT)
            deallocate(UFLXC,_STAT)
            deallocate(DFLXC,_STAT)
            deallocate(DUFLX_DTS,_STAT)
            deallocate(DUFLXC_DTS,_STAT)
            deallocate(CLEARCOUNTS,_STAT)
            deallocate(ALAT,_STAT)
            deallocate(OLRBRG,_STAT)
            deallocate(DOLRBRG_DTS,_STAT)

            call MAPL_TimerOff(MAPL,"--RRTMG",_RC)

         else

            ! Something is wrong. We've selected neither Chou or RRTMG[P]
            _FAIL('No LW radiation code selected!')

         end if SCHEME

         ! Sum up the U and D fluxes to get net downward

         FLX_INT  = FLXD_INT  + FLXU_INT
         FLXA_INT = FLXAD_INT + FLXAU_INT
         FLC_INT  = FLCD_INT  + FLCU_INT
         FLA_INT  = FLAD_INT  + FLAU_INT

         ! Revert to SFCEM to a positive quantity.
         ! Earlier surface emitted positive downwards per Chou-Suarez.
         SFCEM_INT = -SFCEM_INT

         ! RATS Diagnostic <<>> MSL
         if (nRATS .gt. 0) FLX_INT_RAT = FLXD_INT_RAT + FLXU_INT_RAT

         ! Save surface temperature in internal state

         TS_INT    = TS

         ! Export some cloud properties in the infrared

         call MAPL_TimerOn (MAPL,"--MISC")

         call MAPL_GetResource( MAPL, TAUCRIT, 'TAUCRIT:', DEFAULT=0.30, _RC)
         TAUCRIT   = TAUCRIT/2.13

         call MAPL_GetPointer(EXPORT,   CLDPRS,  'CLDPRS'  ,_RC)
         call MAPL_GetPointer(EXPORT,   CLDTMP,  'CLDTMP'  ,_RC)
         call MAPL_GetPointer(EXPORT,    TAUIR,   'TAUIR'  ,_RC)

         if(associated(TAUIR)) TAUIR = 0.5*(TAUDIAG(:,:,:,3)+TAUDIAG(:,:,:,4))

         if(associated(CLDTMP).or.associated(CLDPRS)) then
            if(associated(CLDTMP)) CLDTMP = MAPL_UNDEF
            if(associated(CLDPRS)) CLDPRS = MAPL_UNDEF
            do j=1,jm
               do i=1,im
                  do l=1,lm
                     if(0.5*(TAUDIAG(I,J,L,3)+TAUDIAG(I,J,L,4))>TAUCRIT) then
                        if(associated(CLDTMP)) CLDTMP(I,J) = T  (I,J,L)
                        if(associated(CLDPRS)) CLDPRS(I,J) = PLE(I,J,L-1)
                        exit
                     end if
                  end do
               end do
            end do
         end if

         ! Correcting the timing of the alw and blw (mjs)

         call MAPL_GetPointer(EXPORT,   TSREFF,    'TSREFF' ,_RC)
         call MAPL_GetPointer(EXPORT,   SFCEM,     'SFCEM0' ,_RC)
         call MAPL_GetPointer(EXPORT,   DSFDTS,    'DSFDTS0',_RC)
         call MAPL_GetPointer(EXPORT,   LWS0,      'LWS0'   ,_RC)

         if(associated(TSREFF)) TSREFF = TS             ! reference TS for linearization
         if(associated(DSFDTS)) DSFDTS =-DFDTS(:,:,LM)  ! d(non-negated upward sfc flux) / dTS
         if(associated(SFCEM )) SFCEM  = SFCEM_INT      ! sfc emitted flux (+ve)
         if(associated(LWS0  )) LWS0   = &              ! absorbed (not reflected)
              FLX_INT(:,:,LM) + SFCEM_INT                  !   downward sfc flux (+ve)

         ! Deallocate per-band aerosol arrays

         DEALLOCATE(TAUA)
         DEALLOCATE(SSAA)
         DEALLOCATE(ASYA)

         call MAPL_TimerOff(MAPL,"--MISC")

         !  All done

         RETURN_(ESMF_SUCCESS)

      end subroutine LW_Driver

      ! compute_lw_aer_optics: load and normalize aerosol optical properties
      !   for a column block into aer_props (ty_optical_props_2str).
      !   Called only when need_dirty_optical_props is .true., so aer_props
      !   is guaranteed to be allocated at the call site.
      subroutine compute_lw_aer_optics(colS, colE, &
           TAUA_3d, SSAA_3d, ASYA_3d, aer_props, RC)

         use mo_rte_kind,       only: wp
         use mo_optical_props,  only: ty_optical_props_arry, ty_optical_props_2str

#define TEST_(msg) if (msg /= '') then; write(0,*) trim(msg); VERIFY_(STATUS); end if
         integer,                          intent(in)    :: colS, colE
         real, dimension(:,:,:),           intent(in)    :: TAUA_3d, SSAA_3d, ASYA_3d
         class(ty_optical_props_arry),     intent(inout) :: aer_props
         integer, optional,                intent(out)   :: RC

         integer :: STATUS

         select type (aer_props)
         class is (ty_optical_props_2str)
            ! load unormalized optical properties from aerosol system
            aer_props%tau = real(TAUA_3d(colS:colE,:,:),kind=wp)
            aer_props%ssa = real(SSAA_3d(colS:colE,:,:),kind=wp)
            aer_props%g   = real(ASYA_3d(colS:colE,:,:),kind=wp)
            ! renormalize
            where (aer_props%tau > 0._wp .and. aer_props%ssa > 0._wp )
               aer_props%g   = aer_props%g   / aer_props%ssa
               aer_props%ssa = aer_props%ssa / aer_props%tau
            elsewhere
               aer_props%tau = 0._wp
               aer_props%ssa = 0._wp
               aer_props%g   = 0._wp
            end where

            ! Because RRTMGP is (currently) compiled at R8,
            ! _wp is R8. Apparently with aggressive compiler
            ! flags using Intel, it's possible for, say,
            ! aer_props%ssa to become slightly greater than one
            ! in the above renormalization. So, we add clamps
            ! to the values based on the restrictions see in
            ! RRTMGP/rte-frontend/mo_optical_props.F90
            !
            ! In testing, the values seen were like 1.00000011905028
            ! so just slightly above one.

            ! tau must be greater than 0.0
            aer_props%tau = max(aer_props%tau, 0._wp)
            ! ssa must be between 0.0 and 1.0
            aer_props%ssa = max(min(aer_props%ssa, 1._wp), 0._wp)
            ! g must be between -1.0 and 1.0
            aer_props%g   = max(min(aer_props%g,   1._wp),-1._wp)

         class default
            STATUS = 1
            TEST_('compute_lw_aer_optics: aerosol optical properties hardwired 2-stream for now')
         end select

         RETURN_(ESMF_SUCCESS)
#undef TEST_

      end subroutine compute_lw_aer_optics

      ! compute_lw_cloud_optics_mcica: compute band cloud optical properties,
      !   generate McICA random numbers, sample cloud mask, draw band->gpt,
      !   and apply condensate inhomogeneity scaling.
      subroutine compute_lw_cloud_optics_mcica( &
           colS, colE, ncols_block, LM, ngpt, &
           gen_mro, cond_inhomo, cloud_overlap_type, IM, IM_World, iBeg, jBeg, &
           seeds_time_key, seeds_ctr_key, &
           CWC_3d, REFF_3d, dp_wp, cf_wp, dzmid, &
           cwp_fac_arg, cloud_optics, &
           cloud_props_bnd, cloud_props_gpt, &
           urand, urand_aux, urand_cond, urand_cond_aux, &
           alpha, rcorr, zcw, &
           adl, rdl, &
           cld_mask, &
           MAPL, RC)

         use mo_rte_kind,             only: wp
         use mo_optical_props,        only: ty_optical_props_arry, ty_optical_props_2str
         use mo_cloud_optics_rrtmgp,  only: ty_cloud_optics_rrtmgp
         use mo_cloud_sampling,       only: draw_samples, sampled_mask_max_ran, &
              sampled_urand_gen_max_ran
         use cloud_condensate_inhomogeneity, only: zcw_lookup
#ifdef HAVE_MKL
         use MKL_VSL_TYPE
         use mo_rng_mklvsl_plus,      only: ty_rng_mklvsl_plus
#else
         use mo_rng_mt19937,          only: ty_rng_mt
#endif

#define TEST_(msg) if (msg /= '') then; write(0,*) trim(msg); VERIFY_(STATUS); end if

         integer,                      intent(in)    :: colS, colE, ncols_block, LM, ngpt
         logical,                      intent(in)    :: gen_mro, cond_inhomo
         character(len=*),             intent(in)    :: cloud_overlap_type
         integer,                      intent(in)    :: IM, IM_World, iBeg, jBeg
         integer,                      intent(in)    :: seeds_time_key, seeds_ctr_key
         real, dimension(:,:,:),       intent(in)    :: CWC_3d, REFF_3d
         real(wp), dimension(:,:),     intent(in)    :: dp_wp, cf_wp, dzmid
         real,     dimension(:),       intent(in), optional :: adl, rdl
         real(wp),                     intent(in)    :: cwp_fac_arg
         type(ty_cloud_optics_rrtmgp), intent(inout) :: cloud_optics
         class(ty_optical_props_arry), intent(inout) :: cloud_props_bnd, cloud_props_gpt
         real(wp), dimension(:,:,:),   intent(inout) :: urand
         real(wp), dimension(:,:,:),   intent(inout), optional :: urand_aux
         real(wp), dimension(:,:,:),   intent(inout), optional :: urand_cond, urand_cond_aux
         real(wp), dimension(:,:),     intent(inout), optional :: alpha, rcorr
         real(wp), dimension(:,:,:),   intent(inout), optional :: zcw
         logical,  dimension(:,:,:),   intent(out)   :: cld_mask
         type(MAPL_MetaComp),          intent(inout) :: MAPL
         integer,  optional,           intent(out)   :: RC

         integer :: STATUS
         character(len=256) :: error_msg
         integer :: isub, icol, ilay, igpt, I, J
         integer :: seeds(3)
         real(wp) :: cld_frac
         real :: sigma_qcw
         integer, parameter :: KLIQUID = 2
         integer, parameter :: KICE    = 1
#ifdef HAVE_MKL
         type(ty_rng_mklvsl_plus) :: rng
#else
         type(ty_rng_mt) :: rng
#endif

         ! set PRNG seeds: word1 set per-column below, word2=time, word3=counter
         seeds(2) = seeds_time_key
         seeds(3) = seeds_ctr_key

         !call MAPL_TimerOn(MAPL,"--RRTMGP_CLOUD_OPTICS",RC=STATUS)
         !VERIFY_(STATUS)

         ! Make band in-cloud optical props from cloud_optics and mean in-cloud cloud water paths.
         error_msg = cloud_optics%cloud_optics( &
              real(CWC_3d(colS:colE,:,KLIQUID),kind=wp) * dp_wp(colS:colE,:) * cwp_fac_arg, &
              real(CWC_3d(colS:colE,:,KICE),   kind=wp) * dp_wp(colS:colE,:) * cwp_fac_arg, &
              min( max( real(REFF_3d(colS:colE,:,KLIQUID),kind=wp), &
              cloud_optics%get_min_radius_liq()), cloud_optics%get_max_radius_liq()), &
              min( max( real(REFF_3d(colS:colE,:,KICE),   kind=wp), &
              cloud_optics%get_min_radius_ice()), cloud_optics%get_max_radius_ice()), &
              cloud_props_bnd)
         TEST_(error_msg)

         !call MAPL_TimerOff(MAPL,"--RRTMGP_CLOUD_OPTICS",RC=STATUS)
         !VERIFY_(STATUS)

         !call MAPL_TimerOn(MAPL,"---RRTMGP_MCICA",RC=STATUS)
         !VERIFY_(STATUS)

         ! exponential inter-layer correlations
         if (gen_mro) then
            do ilay = 1,LM-1
               alpha(:,ilay) = exp(-abs(dzmid(colS:colE,ilay))/real(adl(colS:colE),kind=wp))
            enddo
            if (cond_inhomo) then
               do ilay = 1,LM-1
                  rcorr(:,ilay) = exp(-abs(dzmid(colS:colE,ilay))/real(rdl(colS:colE),kind=wp))
               enddo
            endif
         endif

         ! Generate McICA random numbers for block (Philox PRNG)
         do isub = 1, ncols_block
            icol = colS + isub - 1
            J = (icol-1) / IM + 1
            I = icol - (J-1) * IM
            seeds(1) = (jBeg + J - 1) * IM_World + (iBeg + I - 1)
#ifdef HAVE_MKL
            call rng%init(VSL_BRNG_PHILOX4X32X10,seeds)
#else
            call rng%init(seeds)
#endif
            urand(:,:,isub) = reshape(rng%get_random(ngpt*LM),(/ngpt,LM/))
            if (gen_mro) then
               urand_aux(:,:,isub) = reshape(rng%get_random(ngpt*LM),(/ngpt,LM/))
               if (cond_inhomo) then
                  urand_cond    (:,:,isub) = reshape(rng%get_random(ngpt*LM),(/ngpt,LM/))
                  urand_cond_aux(:,:,isub) = reshape(rng%get_random(ngpt*LM),(/ngpt,LM/))
               endif
            endif
            call rng%end()
         end do

         ! cloud sampling to gpoints
         select case (cloud_overlap_type)
         case ("MAX_RAN_OVERLAP")
            error_msg = sampled_mask_max_ran( &
                 urand(:,:,1:ncols_block), cf_wp(colS:colE,:), cld_mask)
            TEST_(error_msg)
         case ("EXP_RAN_OVERLAP")
            STATUS = 1
            TEST_('EXP_RAN_OVERLAP not implemented yet')
         case ("GEN_MAX_RAN_OVERLAP")
            error_msg = sampled_urand_gen_max_ran(alpha, &
                 urand(:,:,1:ncols_block),urand_aux(:,:,1:ncols_block))
            TEST_(error_msg)
            if (cond_inhomo) then
               error_msg = sampled_urand_gen_max_ran(rcorr, &
                    urand_cond(:,:,1:ncols_block),urand_cond_aux(:,:,1:ncols_block))
               TEST_(error_msg)
            end if
            do isub = 1,ncols_block
               icol = colS + isub - 1
               do ilay = 1,LM
                  cld_frac = cf_wp(icol,ilay)
                  if (cld_frac <= 0._wp) then
                     cld_mask(isub,ilay,:) = .false.
                  else
                     cld_mask(isub,ilay,:) = urand(:,ilay,isub) < cld_frac
                     if (cond_inhomo) then
                        if (cld_frac > 0.99_wp) then
                           sigma_qcw = 0.5
                        elseif (cld_frac > 0.9_wp) then
                           sigma_qcw = 0.71
                        else
                           sigma_qcw = 1.0
                        endif
                        do igpt = 1,ngpt
                           if (cld_mask(isub,ilay,igpt)) zcw(isub,ilay,igpt) = &
                                zcw_lookup(real(urand_cond(igpt,ilay,isub)),sigma_qcw)
                        end do
                     end if
                  end if
               end do
            end do
         case default
            STATUS = 1
            TEST_('compute_lw_cloud_optics_mcica: unknown cloud overlap type')
         end select

         ! draw McICA optical property samples (band->gpt)
         TEST_(draw_samples(cld_mask, cloud_props_bnd, cloud_props_gpt))

         ! Apply sub-gridscale condensate scaling
         if (gen_mro) then
            if (cond_inhomo) &
                 where (cld_mask) cloud_props_gpt%tau = cloud_props_gpt%tau * zcw
         end if

         !call MAPL_TimerOff(MAPL,"---RRTMGP_MCICA",RC=STATUS)
         !VERIFY_(STATUS)

         RETURN_(ESMF_SUCCESS)
#undef TEST_

      end subroutine compute_lw_cloud_optics_mcica

      ! compute_lw_gas_optics: compute LW gas optical properties and Planck
      !   source functions for one block of columns.
      subroutine compute_lw_gas_optics(colS, colE, &
           k_dist, p_lay, p_lev, t_lay, t_lev, t_sfc, &
           gas_concs_block, clean_optical_props, sources, &
           MAPL, RC)

         use mo_rte_kind,              only: wp
         use mo_gas_optics_rrtmgp,    only: ty_gas_optics_rrtmgp
         use mo_gas_concentrations,   only: ty_gas_concs
         use mo_optical_props,        only: ty_optical_props_arry
         use mo_source_functions,     only: ty_source_func_lw

#define TEST_(msg) if (msg /= '') then; write(0,*) trim(msg); VERIFY_(STATUS); end if

         integer,                      intent(in)    :: colS, colE
         type(ty_gas_optics_rrtmgp),   intent(inout) :: k_dist
         real(wp), dimension(:,:),     intent(in)    :: p_lay, p_lev, t_lay, t_lev
         real(wp), dimension(:),       intent(in)    :: t_sfc
         type(ty_gas_concs),           intent(inout) :: gas_concs_block
         class(ty_optical_props_arry), intent(inout) :: clean_optical_props
         type(ty_source_func_lw),      intent(inout) :: sources
         type(MAPL_MetaComp),          intent(inout) :: MAPL
         integer, optional,            intent(out)   :: RC

         integer :: STATUS
         character(len=256) :: error_msg

         !call MAPL_TimerOn(MAPL,"---RRTMGP_GAS_OPTICS",RC=STATUS)
         !VERIFY_(STATUS)

         ! get gas optical properties and Planck source functions
         error_msg = k_dist%gas_optics( &
              p_lay(colS:colE,:), p_lev(colS:colE,:), t_lay(colS:colE,:), &
              t_sfc(colS:colE), gas_concs_block, clean_optical_props, sources, &
              tlev = t_lev(colS:colE,:))
         TEST_(error_msg)

         !call MAPL_TimerOff(MAPL,"---RRTMGP_GAS_OPTICS",RC=STATUS)
         !VERIFY_(STATUS)

         RETURN_(ESMF_SUCCESS)
#undef TEST_

      end subroutine compute_lw_gas_optics

      ! compute_lw_rte: solve LW radiative transfer for one block of columns.
      !   Handles clean clear-sky, clean all-sky, dirty clear-sky, and dirty
      !   all-sky cases as controlled by the calc_* / export_* flags.
      subroutine compute_lw_rte( &
           colS, colE, ncols_block, LM, nmom, &
           top_at_1, u2s, nga, &
           calc_clrnoa, calc_allnoa, calc_clrsky, calc_allsky, &
           allnoa_to_allsky_band_xfer_needed, any_band_output, &
           export_clrsky, export_allsky, &
           implements_aerosol_optics, need_dirty_optical_props, &
           clean_optical_props, sources, emis_sfc, &
           dirty_optical_props, aer_props, cloud_props_gpt, &
           flux_up_clrnoa, flux_dn_clrnoa, dfupdts_clrnoa, &
           flux_up_allnoa, flux_dn_allnoa, dfupdts_allnoa, &
           bnd_flux_up_allnoa, bnd_dfupdts_allnoa, &
           flux_up_clrsky, flux_dn_clrsky, dfupdts_clrsky, &
           flux_up_allsky, flux_dn_allsky, dfupdts_allsky, &
           bnd_flux_up_allsky, bnd_dfupdts_allsky, &
           MAPL, RC)

         use mo_rte_kind,        only: wp
         use mo_optical_props,   only: ty_optical_props_arry, ty_optical_props_1scl, &
              ty_optical_props_2str, ty_optical_props_nstr
         use mo_source_functions, only: ty_source_func_lw
         use mo_fluxes,          only: ty_fluxes_broadband
         use mo_fluxes_byband,   only: ty_fluxes_byband
         use mo_rte_lw,          only: rte_lw

#define TEST_(msg) if (msg /= '') then; write(0,*) trim(msg); VERIFY_(STATUS); end if

         integer,                      intent(in)    :: colS, colE, ncols_block, LM, nmom
         logical,                      intent(in)    :: top_at_1, u2s
         integer,                      intent(in)    :: nga
         logical,                      intent(in)    :: calc_clrnoa, calc_allnoa, calc_clrsky, calc_allsky
         logical,                      intent(in)    :: allnoa_to_allsky_band_xfer_needed, any_band_output
         logical,                      intent(in)    :: export_clrsky, export_allsky
         logical,                      intent(in)    :: implements_aerosol_optics, need_dirty_optical_props
         class(ty_optical_props_arry), intent(inout) :: clean_optical_props
         type(ty_source_func_lw),      intent(inout) :: sources
         real(wp), dimension(:,:),     intent(in)    :: emis_sfc
         class(ty_optical_props_arry), intent(inout), optional :: dirty_optical_props
         class(ty_optical_props_arry), intent(inout), optional :: aer_props
         class(ty_optical_props_arry), intent(inout), optional :: cloud_props_gpt
         real(wp), dimension(:,:), intent(inout), target, optional :: flux_up_clrnoa, flux_dn_clrnoa, dfupdts_clrnoa
         real(wp), dimension(:,:), intent(inout), target, optional :: flux_up_allnoa, flux_dn_allnoa, dfupdts_allnoa
         real(wp), dimension(:,:), intent(inout), target, optional :: flux_up_clrsky, flux_dn_clrsky, dfupdts_clrsky
         real(wp), dimension(:,:), intent(inout), target, optional :: flux_up_allsky, flux_dn_allsky, dfupdts_allsky
         real(wp), dimension(:,:,:), intent(inout), target, optional :: bnd_flux_up_allnoa, bnd_dfupdts_allnoa
         real(wp), dimension(:,:,:), intent(inout), target, optional :: bnd_flux_up_allsky, bnd_dfupdts_allsky
         type(MAPL_MetaComp),          intent(inout) :: MAPL
         integer, optional,            intent(out)   :: RC

         integer :: STATUS
         character(len=256) :: error_msg
         type(ty_fluxes_broadband) :: fluxes_clrsky, fluxes_clrnoa, fluxes_allnoa, fluxes_allsky
         type(ty_fluxes_byband)    :: fluxes_byband_allnoa, fluxes_byband_allsky

         !call MAPL_TimerOn(MAPL,"---RRTMGP_RT",RC=STATUS)
         !VERIFY_(STATUS)

         ! clean clear-sky case
         if (calc_clrnoa) then
            fluxes_clrnoa%flux_up     => flux_up_clrnoa(colS:colE,:)
            fluxes_clrnoa%flux_dn     => flux_dn_clrnoa(colS:colE,:)
            fluxes_clrnoa%flux_up_Jac => dfupdts_clrnoa(colS:colE,:)
            error_msg = rte_lw( &
                 clean_optical_props, &
                 top_at_1, sources, emis_sfc(:,colS:colE), &
                 fluxes_clrnoa, n_gauss_angles=nga, use_2stream=u2s)
            TEST_(error_msg)
         end if

         if (present(dirty_optical_props)) then
            ! make copy of clrnoa optical properties as the
            !   starting point for later dirty calculations
            select type (dirty_optical_props)
            class is (ty_optical_props_1scl)
               TEST_(dirty_optical_props%alloc_1scl(ncols_block, LM, clean_optical_props))
            class is (ty_optical_props_2str)
               TEST_(dirty_optical_props%alloc_2str(ncols_block, LM, clean_optical_props))
               select type (clean_optical_props)
               class is (ty_optical_props_2str)
                  dirty_optical_props%ssa = clean_optical_props%ssa
                  dirty_optical_props%g   = clean_optical_props%g
               end select
            class is (ty_optical_props_nstr)
               TEST_(dirty_optical_props%alloc_nstr(nmom, ncols_block, LM, clean_optical_props))
               select type (clean_optical_props)
               class is (ty_optical_props_nstr)
                  dirty_optical_props%ssa = clean_optical_props%ssa
                  dirty_optical_props%p   = clean_optical_props%p
               end select
            end select
            ! all streams have tau
            dirty_optical_props%tau = clean_optical_props%tau
         end if

         ! clean all-sky case
         if (calc_allnoa) then

            ! add in cloud optical properties
            TEST_(cloud_props_gpt%increment(clean_optical_props))

            ! clean all-sky RT
            if (allnoa_to_allsky_band_xfer_needed) then
               fluxes_byband_allnoa%flux_up     => flux_up_allnoa(colS:colE,:)
               fluxes_byband_allnoa%flux_dn     => flux_dn_allnoa(colS:colE,:)
               fluxes_byband_allnoa%flux_up_Jac => dfupdts_allnoa(colS:colE,:)
               fluxes_byband_allnoa%bnd_flux_up     => bnd_flux_up_allnoa(colS:colE,:,:)
               fluxes_byband_allnoa%bnd_flux_up_Jac => bnd_dfupdts_allnoa(colS:colE,:,:)
               error_msg = rte_lw( &
                    clean_optical_props, &
                    top_at_1, sources, emis_sfc(:,colS:colE), &
                    fluxes_byband_allnoa, n_gauss_angles=nga, use_2stream=u2s)
               TEST_(error_msg)
            else
               ! only broadband required
               fluxes_allnoa%flux_up     => flux_up_allnoa(colS:colE,:)
               fluxes_allnoa%flux_dn     => flux_dn_allnoa(colS:colE,:)
               fluxes_allnoa%flux_up_Jac => dfupdts_allnoa(colS:colE,:)
               error_msg = rte_lw( &
                    clean_optical_props, &
                    top_at_1, sources, emis_sfc(:,colS:colE), &
                    fluxes_allnoa, n_gauss_angles=nga, use_2stream=u2s)
               TEST_(error_msg)
            endif
         end if

         if (export_clrsky .or. export_allsky) then
            if (implements_aerosol_optics) then

               ! dirty flux calculations required ...

               ! "dirty_optical_props" is currently just a copy of the clrnoa optical_props
               !   so must now add in aerosols to make it actually dirty
               TEST_(aer_props%increment(dirty_optical_props))

               ! dirty clear-sky RT
               if (calc_clrsky) then
                  fluxes_clrsky%flux_up     => flux_up_clrsky(colS:colE,:)
                  fluxes_clrsky%flux_dn     => flux_dn_clrsky(colS:colE,:)
                  fluxes_clrsky%flux_up_Jac => dfupdts_clrsky(colS:colE,:)
                  error_msg = rte_lw( &
                       dirty_optical_props, &
                       top_at_1, sources, emis_sfc(:,colS:colE), &
                       fluxes_clrsky, n_gauss_angles=nga, use_2stream=u2s)
                  TEST_(error_msg)
               end if

               ! dirty all-sky case
               if (calc_allsky) then

                  ! add in cloud optical properties
                  TEST_(cloud_props_gpt%increment(dirty_optical_props))

                  ! dirty all-sky RT
                  ! (band output currently only available for all-sky case)
                  if (any_band_output) then
                     fluxes_byband_allsky%flux_up     => flux_up_allsky(colS:colE,:)
                     fluxes_byband_allsky%flux_dn     => flux_dn_allsky(colS:colE,:)
                     fluxes_byband_allsky%flux_up_Jac => dfupdts_allsky(colS:colE,:)
                     fluxes_byband_allsky%bnd_flux_up     => bnd_flux_up_allsky(colS:colE,:,:)
                     fluxes_byband_allsky%bnd_flux_up_Jac => bnd_dfupdts_allsky(colS:colE,:,:)
                     error_msg = rte_lw( &
                          dirty_optical_props, &
                          top_at_1, sources, emis_sfc(:,colS:colE), &
                          fluxes_byband_allsky, n_gauss_angles=nga, use_2stream=u2s)
                     TEST_(error_msg)
                  else
                     fluxes_allsky%flux_up     => flux_up_allsky(colS:colE,:)
                     fluxes_allsky%flux_dn     => flux_dn_allsky(colS:colE,:)
                     fluxes_allsky%flux_up_Jac => dfupdts_allsky(colS:colE,:)
                     error_msg = rte_lw( &
                          dirty_optical_props, &
                          top_at_1, sources, emis_sfc(:,colS:colE), &
                          fluxes_allsky, n_gauss_angles=nga, use_2stream=u2s)
                     TEST_(error_msg)
                  end if
               end if

            else

               ! there are no aerosols so we are done because the
               !   dirty cases are the same as the clean ones
               if (export_clrsky) then
                  flux_up_clrsky(colS:colE,:) = flux_up_clrnoa(colS:colE,:)
                  flux_dn_clrsky(colS:colE,:) = flux_dn_clrnoa(colS:colE,:)
                  dfupdts_clrsky(colS:colE,:) = dfupdts_clrnoa(colS:colE,:)
               end if
               if (export_allsky) then
                  flux_up_allsky(colS:colE,:) = flux_up_allnoa(colS:colE,:)
                  flux_dn_allsky(colS:colE,:) = flux_dn_allnoa(colS:colE,:)
                  dfupdts_allsky(colS:colE,:) = dfupdts_allnoa(colS:colE,:)
                  if (any_band_output) then
                     bnd_flux_up_allsky(colS:colE,:,:) = bnd_flux_up_allnoa(colS:colE,:,:)
                     bnd_dfupdts_allsky(colS:colE,:,:) = bnd_dfupdts_allnoa(colS:colE,:,:)
                  end if
               end if

            end if ! implements_aerosol_optics
         end if ! export dirty clear-sky or all-sky

         !call MAPL_TimerOff(MAPL,"---RRTMGP_RT",RC=STATUS)
         !VERIFY_(STATUS)

         RETURN_(ESMF_SUCCESS)
#undef TEST_

      end subroutine compute_lw_rte

      ! PROCESS_RRTMGP_LW_BLOCK: process one block of columns through the
      !   full LW RRTMGP pipeline (aerosol optics, cloud optics, gas optics,
      !   RTE solve).  Intended to be called from a serial or OpenMP
      !   parallel do loop over blocks.
      subroutine PROCESS_RRTMGP_LW_BLOCK( &
           b, rrtmgp_blockSize, ncol, LM, nmom, ngpt, nga, &
           IM, IM_World, iBeg, jBeg, &
           top_at_1, u2s, &
           seeds_time_key, seeds_ctr_key, &
           cwp_fac, &
           need_cloud_optical_props, need_dirty_optical_props, &
           gen_mro, cond_inhomo, cloud_overlap_type, &
           calc_clrnoa, calc_allnoa, calc_clrsky, calc_allsky, &
           allnoa_to_allsky_band_xfer_needed, any_band_output, &
           export_clrsky, export_allsky, implements_aerosol_optics, &
           k_dist, cloud_optics, gas_concs, &
           p_lay, p_lev, t_lay, t_lev, t_sfc, dp_wp, cf_wp, dzmid, emis_sfc, &
           adl, rdl, &
           CWC_3d, REFF_3d, &
           TAUA_3d, SSAA_3d, ASYA_3d, &
           flux_up_clrnoa, flux_dn_clrnoa, dfupdts_clrnoa, &
           flux_up_allnoa, flux_dn_allnoa, dfupdts_allnoa, &
           bnd_flux_up_allnoa, bnd_dfupdts_allnoa, &
           flux_up_clrsky, flux_dn_clrsky, dfupdts_clrsky, &
           flux_up_allsky, flux_dn_allsky, dfupdts_allsky, &
           bnd_flux_up_allsky, bnd_dfupdts_allsky, &
           MAPL, RC)

         use mo_rte_kind,             only: wp
         use mo_gas_optics_rrtmgp,   only: ty_gas_optics_rrtmgp
         use mo_gas_concentrations,  only: ty_gas_concs
         use mo_optical_props,       only: ty_optical_props_2str
         use mo_source_functions,    only: ty_source_func_lw
         use mo_cloud_optics_rrtmgp, only: ty_cloud_optics_rrtmgp

#define TEST_(msg) if (msg /= '') then; write(0,*) trim(msg); VERIFY_(STATUS); end if

         integer,                      intent(in)    :: b, rrtmgp_blockSize, ncol, LM, nmom, ngpt, nga
         integer,                      intent(in)    :: IM, IM_World, iBeg, jBeg
         logical,                      intent(in)    :: top_at_1, u2s
         integer,                      intent(in)    :: seeds_time_key, seeds_ctr_key
         real(wp),                     intent(in)    :: cwp_fac
         logical,                      intent(in)    :: need_cloud_optical_props, need_dirty_optical_props
         logical,                      intent(in)    :: gen_mro, cond_inhomo
         character(len=*),             intent(in)    :: cloud_overlap_type
         logical,                      intent(in)    :: calc_clrnoa, calc_allnoa, calc_clrsky, calc_allsky
         logical,                      intent(in)    :: allnoa_to_allsky_band_xfer_needed, any_band_output
         logical,                      intent(in)    :: export_clrsky, export_allsky, implements_aerosol_optics
         type(ty_gas_optics_rrtmgp),   intent(inout) :: k_dist
         type(ty_cloud_optics_rrtmgp), intent(inout) :: cloud_optics
         type(ty_gas_concs),           intent(inout) :: gas_concs
         real(wp), dimension(:,:),     intent(in)    :: p_lay, p_lev, t_lay, t_lev
         real(wp), dimension(:),       intent(in)    :: t_sfc
         real(wp), dimension(:,:),     intent(in)    :: dp_wp, cf_wp, dzmid
         real(wp), dimension(:,:),     intent(in)    :: emis_sfc
         real,     dimension(:),       intent(in), optional :: adl, rdl
         real,     dimension(:,:,:),   pointer       :: CWC_3d, REFF_3d
         real,     dimension(:,:,:),   pointer       :: TAUA_3d, SSAA_3d, ASYA_3d
         real(wp), dimension(:,:), intent(inout), target, optional :: flux_up_clrnoa, flux_dn_clrnoa, dfupdts_clrnoa
         real(wp), dimension(:,:), intent(inout), target, optional :: flux_up_allnoa, flux_dn_allnoa, dfupdts_allnoa
         real(wp), dimension(:,:), intent(inout), target, optional :: flux_up_clrsky, flux_dn_clrsky, dfupdts_clrsky
         real(wp), dimension(:,:), intent(inout), target, optional :: flux_up_allsky, flux_dn_allsky, dfupdts_allsky
         real(wp), dimension(:,:,:), intent(inout), target, optional :: bnd_flux_up_allnoa, bnd_dfupdts_allnoa
         real(wp), dimension(:,:,:), intent(inout), target, optional :: bnd_flux_up_allsky, bnd_dfupdts_allsky
         type(MAPL_MetaComp),          intent(inout) :: MAPL
         integer,  optional,           intent(out)   :: RC

         integer :: STATUS
         character(len=256) :: error_msg
         integer :: ncols_block, colS, colE

         ! local RRTMGP objects (LW always uses 2-stream)
         type(ty_optical_props_2str) :: clean_optical_props
         type(ty_optical_props_2str) :: dirty_optical_props
         type(ty_optical_props_2str) :: aer_props
         type(ty_optical_props_2str) :: cloud_props_bnd, cloud_props_gpt
         type(ty_source_func_lw)     :: sources
         type(ty_gas_concs)          :: gas_concs_block

         ! per-block scratch arrays
         real(wp), dimension(:,:,:), allocatable :: urand, urand_aux, urand_cond, urand_cond_aux
         real(wp), dimension(:,:,:), allocatable :: zcw
         real(wp), dimension(:,:),   allocatable :: alpha, rcorr
         logical,  dimension(:,:,:), allocatable :: cld_mask

         ! compute column range for this block (final block may be partial)
         ncols_block = min(rrtmgp_blockSize, ncol - (b-1)*rrtmgp_blockSize)
         colS = (b-1) * rrtmgp_blockSize + 1
         colE = colS + ncols_block - 1

         ! spectral init + array allocation for gas optics and Planck sources
         TEST_(clean_optical_props%init(k_dist))
         TEST_(clean_optical_props%alloc_2str(ncols_block, LM))
         TEST_(sources%init(k_dist))
         TEST_(sources%alloc(ncols_block, LM))

         ! subset gas concentrations for this block
         TEST_(gas_concs%get_subset(colS, ncols_block, gas_concs_block))

         ! aerosol optics objects (always 2-stream for LW)
         if (need_dirty_optical_props) then
            TEST_(dirty_optical_props%init(k_dist))
            TEST_(aer_props%init(k_dist%get_band_lims_wavenumber()))
            TEST_(aer_props%alloc_2str(ncols_block, LM))
         end if

         ! cloud optics objects and scratch arrays
         if (need_cloud_optical_props) then
            TEST_(cloud_props_bnd%init(k_dist%get_band_lims_wavenumber()))
            TEST_(cloud_props_bnd%alloc_2str(ncols_block, LM))
            TEST_(cloud_props_gpt%init(k_dist))
            TEST_(cloud_props_gpt%alloc_2str(ncols_block, LM))
            allocate(urand(ngpt, LM, ncols_block), _STAT)
            allocate(cld_mask(ncols_block, LM, ngpt), _STAT)
            if (gen_mro) then
               allocate(urand_aux(ngpt, LM, ncols_block), _STAT)
               allocate(alpha(ncols_block, LM-1), _STAT)
               if (cond_inhomo) then
                  allocate(urand_cond    (ngpt, LM, ncols_block), _STAT)
                  allocate(urand_cond_aux(ngpt, LM, ncols_block), _STAT)
                  allocate(rcorr(ncols_block, LM-1), _STAT)
                  allocate(zcw  (ncols_block, LM, ngpt), _STAT)
               end if
            end if
         end if

         ! aerosol optical properties
         if (need_dirty_optical_props) then
            call compute_lw_aer_optics(colS, colE, &
                 TAUA_3d, SSAA_3d, ASYA_3d, aer_props, _RC)
         end if

         ! cloud optical properties (McICA sampling)
         if (need_cloud_optical_props) then
            call compute_lw_cloud_optics_mcica( &
                 colS, colE, ncols_block, LM, ngpt, &
                 gen_mro, cond_inhomo, cloud_overlap_type, IM, IM_World, iBeg, jBeg, &
                 seeds_time_key, seeds_ctr_key, &
                 CWC_3d, REFF_3d, dp_wp, cf_wp, dzmid, &
                 cwp_fac, cloud_optics, &
                 cloud_props_bnd, cloud_props_gpt, &
                 urand, &
                 urand_aux=urand_aux, urand_cond=urand_cond, urand_cond_aux=urand_cond_aux, &
                 alpha=alpha, rcorr=rcorr, zcw=zcw, &
                 adl=adl, rdl=rdl, &
                 cld_mask=cld_mask, &
                 MAPL=MAPL, _RC)
         end if

         ! gas optical properties and Planck source functions
         call compute_lw_gas_optics(colS, colE, &
              k_dist, p_lay, p_lev, t_lay, t_lev, t_sfc, &
              gas_concs_block, clean_optical_props, sources, &
              MAPL=MAPL, _RC)

         ! radiative transfer solve (conditional on which optional objects are present)
         if (need_dirty_optical_props .and. need_cloud_optical_props) then
            call compute_lw_rte( &
                 colS, colE, ncols_block, LM, nmom, &
                 top_at_1, u2s, nga, &
                 calc_clrnoa, calc_allnoa, calc_clrsky, calc_allsky, &
                 allnoa_to_allsky_band_xfer_needed, any_band_output, &
                 export_clrsky, export_allsky, &
                 implements_aerosol_optics, need_dirty_optical_props, &
                 clean_optical_props, sources, emis_sfc, &
                 dirty_optical_props=dirty_optical_props, aer_props=aer_props, &
                 cloud_props_gpt=cloud_props_gpt, &
                 flux_up_clrnoa=flux_up_clrnoa, flux_dn_clrnoa=flux_dn_clrnoa, dfupdts_clrnoa=dfupdts_clrnoa, &
                 flux_up_allnoa=flux_up_allnoa, flux_dn_allnoa=flux_dn_allnoa, dfupdts_allnoa=dfupdts_allnoa, &
                 bnd_flux_up_allnoa=bnd_flux_up_allnoa, bnd_dfupdts_allnoa=bnd_dfupdts_allnoa, &
                 flux_up_clrsky=flux_up_clrsky, flux_dn_clrsky=flux_dn_clrsky, dfupdts_clrsky=dfupdts_clrsky, &
                 flux_up_allsky=flux_up_allsky, flux_dn_allsky=flux_dn_allsky, dfupdts_allsky=dfupdts_allsky, &
                 bnd_flux_up_allsky=bnd_flux_up_allsky, bnd_dfupdts_allsky=bnd_dfupdts_allsky, &
                 MAPL=MAPL, RC=STATUS)
         else if (need_dirty_optical_props) then
            call compute_lw_rte( &
                 colS, colE, ncols_block, LM, nmom, &
                 top_at_1, u2s, nga, &
                 calc_clrnoa, calc_allnoa, calc_clrsky, calc_allsky, &
                 allnoa_to_allsky_band_xfer_needed, any_band_output, &
                 export_clrsky, export_allsky, &
                 implements_aerosol_optics, need_dirty_optical_props, &
                 clean_optical_props, sources, emis_sfc, &
                 dirty_optical_props=dirty_optical_props, aer_props=aer_props, &
                 flux_up_clrnoa=flux_up_clrnoa, flux_dn_clrnoa=flux_dn_clrnoa, dfupdts_clrnoa=dfupdts_clrnoa, &
                 flux_up_allnoa=flux_up_allnoa, flux_dn_allnoa=flux_dn_allnoa, dfupdts_allnoa=dfupdts_allnoa, &
                 bnd_flux_up_allnoa=bnd_flux_up_allnoa, bnd_dfupdts_allnoa=bnd_dfupdts_allnoa, &
                 flux_up_clrsky=flux_up_clrsky, flux_dn_clrsky=flux_dn_clrsky, dfupdts_clrsky=dfupdts_clrsky, &
                 flux_up_allsky=flux_up_allsky, flux_dn_allsky=flux_dn_allsky, dfupdts_allsky=dfupdts_allsky, &
                 bnd_flux_up_allsky=bnd_flux_up_allsky, bnd_dfupdts_allsky=bnd_dfupdts_allsky, &
                 MAPL=MAPL, RC=STATUS)
         else if (need_cloud_optical_props) then
            call compute_lw_rte( &
                 colS, colE, ncols_block, LM, nmom, &
                 top_at_1, u2s, nga, &
                 calc_clrnoa, calc_allnoa, calc_clrsky, calc_allsky, &
                 allnoa_to_allsky_band_xfer_needed, any_band_output, &
                 export_clrsky, export_allsky, &
                 implements_aerosol_optics, need_dirty_optical_props, &
                 clean_optical_props, sources, emis_sfc, &
                 cloud_props_gpt=cloud_props_gpt, &
                 flux_up_clrnoa=flux_up_clrnoa, flux_dn_clrnoa=flux_dn_clrnoa, dfupdts_clrnoa=dfupdts_clrnoa, &
                 flux_up_allnoa=flux_up_allnoa, flux_dn_allnoa=flux_dn_allnoa, dfupdts_allnoa=dfupdts_allnoa, &
                 bnd_flux_up_allnoa=bnd_flux_up_allnoa, bnd_dfupdts_allnoa=bnd_dfupdts_allnoa, &
                 flux_up_clrsky=flux_up_clrsky, flux_dn_clrsky=flux_dn_clrsky, dfupdts_clrsky=dfupdts_clrsky, &
                 flux_up_allsky=flux_up_allsky, flux_dn_allsky=flux_dn_allsky, dfupdts_allsky=dfupdts_allsky, &
                 bnd_flux_up_allsky=bnd_flux_up_allsky, bnd_dfupdts_allsky=bnd_dfupdts_allsky, &
                 MAPL=MAPL, RC=STATUS)
         else
            call compute_lw_rte( &
                 colS, colE, ncols_block, LM, nmom, &
                 top_at_1, u2s, nga, &
                 calc_clrnoa, calc_allnoa, calc_clrsky, calc_allsky, &
                 allnoa_to_allsky_band_xfer_needed, any_band_output, &
                 export_clrsky, export_allsky, &
                 implements_aerosol_optics, need_dirty_optical_props, &
                 clean_optical_props, sources, emis_sfc, &
                 flux_up_clrnoa=flux_up_clrnoa, flux_dn_clrnoa=flux_dn_clrnoa, dfupdts_clrnoa=dfupdts_clrnoa, &
                 flux_up_allnoa=flux_up_allnoa, flux_dn_allnoa=flux_dn_allnoa, dfupdts_allnoa=dfupdts_allnoa, &
                 bnd_flux_up_allnoa=bnd_flux_up_allnoa, bnd_dfupdts_allnoa=bnd_dfupdts_allnoa, &
                 flux_up_clrsky=flux_up_clrsky, flux_dn_clrsky=flux_dn_clrsky, dfupdts_clrsky=dfupdts_clrsky, &
                 flux_up_allsky=flux_up_allsky, flux_dn_allsky=flux_dn_allsky, dfupdts_allsky=dfupdts_allsky, &
                 bnd_flux_up_allsky=bnd_flux_up_allsky, bnd_dfupdts_allsky=bnd_dfupdts_allsky, &
                 MAPL=MAPL, RC=STATUS)
         end if
         VERIFY_(STATUS)

         ! finalize/deallocate per-block RRTMGP objects
         call sources%finalize()
         call clean_optical_props%finalize()
         if (need_dirty_optical_props) then
            call dirty_optical_props%finalize()
            call aer_props%finalize()
         end if
         if (need_cloud_optical_props) then
            call cloud_props_bnd%finalize()
            call cloud_props_gpt%finalize()
            deallocate(urand, cld_mask, _STAT)
            if (gen_mro) then
               deallocate(urand_aux, alpha, _STAT)
               if (cond_inhomo) then
                  deallocate(urand_cond, urand_cond_aux, rcorr, zcw, _STAT)
               end if
            end if
         end if

         RETURN_(ESMF_SUCCESS)
#undef TEST_

      end subroutine PROCESS_RRTMGP_LW_BLOCK

      subroutine Update_Flx(IM,JM,LM,RC)
         use mo_rte_kind, only: wp
         integer,           intent(IN ) :: IM, JM, LM
         integer, optional, intent(OUT) :: RC

         !  Locals

         character(len=ESMF_MAXSTR)        :: Iam
         integer                           :: STATUS

         real,          dimension(IM,JM)   :: DELT
         integer                           :: K
         integer                           :: N !<<>> MSL
         integer                           :: LEV_LOW_MID
         integer                           :: LEV_MID_HIGH
         real                              :: PRS_LOW_MID                   ! pressure separating low and middle clouds
         real                              :: PRS_MID_HIGH                  ! pressure separating low and high   clouds

         ! band wavenumber bounds (m-1)
         real :: wn1, wn2

         ! pointer to import

         real, pointer, dimension(:,:  )   :: TSINST

         ! pointers to export

         real, pointer, dimension(:,:,:)   :: FLX
         real, pointer, dimension(:,:,:)   :: FLXA
         real, pointer, dimension(:,:,:)   :: FLC
         real, pointer, dimension(:,:,:)   :: FLA
         real, pointer, dimension(:,:,:)   :: FLXU
         real, pointer, dimension(:,:,:)   :: FLXAU
         real, pointer, dimension(:,:,:)   :: FLCU
         real, pointer, dimension(:,:,:)   :: FLAU
         real, pointer, dimension(:,:,:)   :: FLXD
         real, pointer, dimension(:,:,:)   :: FLXAD
         real, pointer, dimension(:,:,:)   :: FLCD
         real, pointer, dimension(:,:,:)   :: FLAD
         real, pointer, dimension(:,:  )   :: TSREFF
         real, pointer, dimension(:,:  )   :: SFCEM
         real, pointer, dimension(:,:  )   :: DSFDTS
         real, pointer, dimension(:,:  )   :: SFCEM0
         real, pointer, dimension(:,:  )   :: DSFDTS0
         real, pointer, dimension(:,:  )   :: OLR
         real, pointer, dimension(:,:  )   :: OLRA
         real, pointer, dimension(:,:  )   :: OLC
         real, pointer, dimension(:,:  )   :: OLCC5
         real, pointer, dimension(:,:  )   :: OLA
         real, pointer, dimension(:,:  )   :: FLNS
         real, pointer, dimension(:,:  )   :: FLNSNA
         real, pointer, dimension(:,:  )   :: FLNSC
         real, pointer, dimension(:,:  )   :: FLNSA
         real, pointer, dimension(:,:  )   :: LWS
         real, pointer, dimension(:,:  )   :: LWSA
         real, pointer, dimension(:,:  )   :: LCS
         real, pointer, dimension(:,:  )   :: LCSC5
         real, pointer, dimension(:,:  )   :: LAS
         real, pointer, dimension(:,:  )   :: CLDTT
         real, pointer, dimension(:,:  )   :: ptr2d

         real, pointer, dimension(:,:,:)   :: FCLD
         real, pointer, dimension(:    )   :: PREF

         real, allocatable, dimension(:,:) :: DUMTT, OLRB

         ! RATS diagnostics <<>> MSL
         real, pointer, dimension(:,:  )   :: RAT_2D, EMIS
         real, pointer, dimension(:,:,:)   :: RAT_3D

         ! access to RRTMGP wavenumber limits
         real(wp) :: band_lims_wvn(2,nbndlw)

         !  Begin...

         IAm = "Update_Flx"

         ! Pointers to Exports

         call MAPL_GetPointer(EXPORT,   FLX   ,    'FLX',   _RC)
         call MAPL_GetPointer(EXPORT,   FLXA  ,    'FLXA',  _RC)
         call MAPL_GetPointer(EXPORT,   FLC   ,    'FLC',   _RC)
         call MAPL_GetPointer(EXPORT,   FLA   ,    'FLA',   _RC)
         call MAPL_GetPointer(EXPORT,   FLXU  ,    'FLXU',  _RC)
         call MAPL_GetPointer(EXPORT,   FLXAU ,    'FLXAU', _RC)
         call MAPL_GetPointer(EXPORT,   FLCU  ,    'FLCU',  _RC)
         call MAPL_GetPointer(EXPORT,   FLAU  ,    'FLAU',  _RC)
         call MAPL_GetPointer(EXPORT,   FLXD  ,    'FLXD',  _RC)
         call MAPL_GetPointer(EXPORT,   FLXAD ,    'FLXAD', _RC)
         call MAPL_GetPointer(EXPORT,   FLCD  ,    'FLCD',  _RC)
         call MAPL_GetPointer(EXPORT,   FLAD  ,    'FLAD',  _RC)
         call MAPL_GetPointer(EXPORT,   TSREFF,    'TSREFF',_RC)
         call MAPL_GetPointer(EXPORT,   SFCEM ,    'SFCEM', _RC)
         call MAPL_GetPointer(EXPORT,   DSFDTS,    'DSFDTS',_RC)
         call MAPL_GetPointer(EXPORT,   SFCEM0,    'SFCEM0',_RC)
         call MAPL_GetPointer(EXPORT,  DSFDTS0,   'DSFDTS0',_RC)
         call MAPL_GetPointer(EXPORT,   OLR   ,    'OLR'   ,_RC)
         call MAPL_GetPointer(EXPORT,   OLRA  ,    'OLRA'  ,_RC)
         call MAPL_GetPointer(EXPORT,   OLC   ,    'OLC'   ,_RC)
         call MAPL_GetPointer(EXPORT,   OLCC5 ,    'OLCC5' ,_RC)
         call MAPL_GetPointer(EXPORT,   OLA   ,    'OLA'   ,_RC)
         call MAPL_GetPointer(EXPORT,   LWS   ,    'LWS'   ,_RC)
         call MAPL_GetPointer(EXPORT,   LWSA  ,    'LWSA'  ,_RC)
         call MAPL_GetPointer(EXPORT,   LCS   ,    'LCS'   ,_RC)
         call MAPL_GetPointer(EXPORT,   LCSC5 ,    'LCSC5' ,_RC)
         call MAPL_GetPointer(EXPORT,   LAS   ,    'LAS'   ,_RC)
         call MAPL_GetPointer(EXPORT,   FLNS  ,    'FLNS'  ,_RC)
         call MAPL_GetPointer(EXPORT,   FLNSNA,    'FLNSNA',_RC)
         call MAPL_GetPointer(EXPORT,   FLNSC ,    'FLNSC' ,_RC)
         call MAPL_GetPointer(EXPORT,   FLNSA ,    'FLNSA' ,_RC)

         call MAPL_GetPointer(EXPORT,   CLDTT ,  'CLDTT'   ,ALLOC=.TRUE.,_RC)

         ! Determine the 2-D Total Cloud Fraction

         call MAPL_GetResource( MAPL, PRS_LOW_MID,    'PRS_LOW_MID_CLOUDS:' ,   DEFAULT=70000.,      _RC)
         call MAPL_GetResource( MAPL, PRS_MID_HIGH,   'PRS_MID_HIGH_CLOUDS:',   DEFAULT=40000.,      _RC)

         call MAPL_GetPointer( IMPORT, FCLD, 'FCLD', _RC)
         call MAPL_GetPointer( IMPORT, PREF, 'PREF', _RC)

         ALLOCATE( DUMTT(IM,JM), STAT=STATUS)
         VERIFY_(STATUS)

         ! Determine the model level separating mid and high clouds
         LEV_MID_HIGH = 1
         do K = 1, LM
            if( PREF(K) >= PRS_MID_HIGH ) then
               LEV_MID_HIGH = K
               exit
            end if
         end do

         ! Determine the model level seperating low and middle clouds
         LEV_LOW_MID = LM
         do K = 1, LM
            if( PREF(K) >= PRS_LOW_MID  ) then
               LEV_LOW_MID = K
               exit
            end if
         end do

         DUMTT = 0.
         do K=1,LEV_MID_HIGH-1
            DUMTT = max(DUMTT,FCLD(:,:,K))
         end do
         CLDTT = (1-DUMTT)
         DUMTT = 0.
         do K= LEV_MID_HIGH,LEV_LOW_MID-1
            DUMTT = max(DUMTT,FCLD(:,:,K))
         end do
         CLDTT = CLDTT*(1-DUMTT)
         DUMTT = 0.
         do K=LEV_LOW_MID,LM
            DUMTT = max(DUMTT,FCLD(:,:,K))
         end do
         CLDTT = 1.0 - CLDTT*(1-DUMTT)

         ! Pointers to Imports

         call MAPL_GetPointer(IMPORT,   TSINST, 'TSINST',   _RC)

         ! Update fluxes

         ! linearization with surface temperature notes:
         ! a. only upward fluxes linearized wrt surface temperature
         ! b. the derivatives DFDTS[C] have the same sign convention as the negated upward fluxes
         !      (i.e., they are the derivatives of negated upward fluxes with surface temperature)

         ! surface temperature change since refresh for linearization
         DELT = TSINST - TS_INT

         if( USE_CHOU .or. USE_RRTMGP ) THEN

            ! fill 3D fluxes
            do K = 0, LM
               ! net downward (downward plus negated upward) fluxes
               if(associated(FLX))    FLX (:,:,K) =   FLX_INT(:,:,K) + DFDTS   (:,:,K) * DELT ! all-sky
               if(associated(FLXA))  FLXA (:,:,K) =  FLXA_INT(:,:,K) + DFDTSNA (:,:,K) * DELT ! all-sky no-aerosol
               if(associated(FLC))    FLC (:,:,K) =   FLC_INT(:,:,K) + DFDTSC  (:,:,K) * DELT ! clr-sky
               if(associated(FLA))    FLA (:,:,K) =   FLA_INT(:,:,K) + DFDTSCNA(:,:,K) * DELT ! clr-sky no-aerosol
               ! negated upward fluxes
               if(associated(FLXU))   FLXU(:,:,K) =  FLXU_INT(:,:,K) + DFDTS   (:,:,K) * DELT
               if(associated(FLXAU)) FLXAU(:,:,K) = FLXAU_INT(:,:,K) + DFDTSNA (:,:,K) * DELT
               if(associated(FLCU))   FLCU(:,:,K) =  FLCU_INT(:,:,K) + DFDTSC  (:,:,K) * DELT
               if(associated(FLAU))   FLAU(:,:,K) =  FLAU_INT(:,:,K) + DFDTSCNA(:,:,K) * DELT
               ! downward fluxes
               if(associated(FLXD))   FLXD(:,:,K) =  FLXD_INT(:,:,K)
               if(associated(FLXAD)) FLXAD(:,:,K) = FLXAD_INT(:,:,K)
               if(associated(FLCD))   FLCD(:,:,K) =  FLCD_INT(:,:,K)
               if(associated(FLAD))   FLAD(:,:,K) =  FLAD_INT(:,:,K)
            end do

            ! fill TOA exports
            ! outgoing longwave radiation
            ! pmn: using FLXU_INT, etc. would be better ... here assuming down at TOA is zero
            if(associated(OLR  )) OLR   = -( FLX_INT(:,:, 0) + DFDTS   (:,:, 0) * DELT)
            if(associated(OLRA )) OLRA  = -(FLXA_INT(:,:, 0) + DFDTSNA (:,:, 0) * DELT)
            if(associated(OLC  )) OLC   = -( FLC_INT(:,:, 0) + DFDTSC  (:,:, 0) * DELT)
            if(associated(OLA  )) OLA   = -( FLA_INT(:,:, 0) + DFDTSCNA(:,:, 0) * DELT)
            if(associated(OLCC5)) then
               where(CLDTT <= 0.05 )
                  OLCC5 = -( FLC_INT(:,:, 0) + DFDTSC  (:,:, 0) * DELT)
               elsewhere
                  OLCC5 = MAPL_UNDEF
               endwhere
            endif

            ! fill surface exports

            ! current surface emitted flux derivative wrt surface temperature (+ve)
            ! pmn: should be deprecated ... same as DSFDTS0
            if(associated(DSFDTS)) DSFDTS = -DFDTS(:,:,LM)

            ! surface emitted flux (+ve)
            if(associated(SFCEM)) SFCEM = SFCEM_INT - DFDTS(:,:,LM) * DELT

            ! absorbed (non-reflected) downward surface fluxes
            ! (remember: downward fluxes are not not linearized)
            if(associated(LWS  )) LWS   =  FLX_INT(:,:,LM) + SFCEM_INT
            if(associated(LWSA )) LWSA  = FLXA_INT(:,:,LM) + SFCEM_INT
            if(associated(LCS  )) LCS   =  FLC_INT(:,:,LM) + SFCEM_INT
            if(associated(LAS  )) LAS   =  FLA_INT(:,:,LM) + SFCEM_INT
            if(associated(LCSC5)) then
               where(CLDTT <= 0.05 )
                  LCSC5 =  FLC_INT(:,:,LM) + SFCEM_INT
               elsewhere
                  LCSC5 = MAPL_UNDEF
               endwhere
            endif

            ! surface net downward fluxes
            if(associated(FLNS  )) FLNS   =  FLX_INT(:,:,LM) + DFDTS   (:,:,LM) * DELT
            if(associated(FLNSNA)) FLNSNA = FLXA_INT(:,:,LM) + DFDTSNA (:,:,LM) * DELT
            if(associated(FLNSC )) FLNSC  =  FLC_INT(:,:,LM) + DFDTSC  (:,:,LM) * DELT
            if(associated(FLNSA )) FLNSA  =  FLA_INT(:,:,LM) + DFDTSCNA(:,:,LM) * DELT

            ! RRTMG is a special case because its no-aerosol cases are missing
         else if( USE_RRTMG ) THEN

            ! fill 3D fluxes
            do K = 0, LM
               ! net downward (downward plus negated upward) fluxes
               if(associated(FLX))    FLX (:,:,K) =  FLX_INT(:,:,K) + DFDTS (:,:,K) * DELT ! all-sky
               if(associated(FLXA))  FLXA (:,:,K) = MAPL_UNDEF                             ! all-sky no-aerosol
               if(associated(FLC))    FLC (:,:,K) =  FLC_INT(:,:,K) + DFDTSC(:,:,K) * DELT ! clr-sky
               if(associated(FLA))    FLA (:,:,K) = MAPL_UNDEF                             ! clr-sky no-aerosol
               ! negated upward fluxes
               if(associated(FLXU))   FLXU(:,:,K) = FLXU_INT(:,:,K) + DFDTS (:,:,K) * DELT
               if(associated(FLXAU)) FLXAU(:,:,K) = MAPL_UNDEF
               if(associated(FLCU))   FLCU(:,:,K) = FLCU_INT(:,:,K) + DFDTSC(:,:,K) * DELT
               if(associated(FLAU))   FLAU(:,:,K) = MAPL_UNDEF
               ! downward fluxes
               if(associated(FLXD))   FLXD(:,:,K) = FLXD_INT(:,:,K)
               if(associated(FLXAD)) FLXAD(:,:,K) = MAPL_UNDEF
               if(associated(FLCD))   FLCD(:,:,K) = FLCD_INT(:,:,K)
               if(associated(FLAD))   FLAD(:,:,K) = MAPL_UNDEF
            end do

            ! fill TOA exports
            ! outgoing longwave radiation
            ! pmn: using FLXU_INT, etc. would be better ... here assuming down at TOA is zero
            if(associated(OLR  )) OLR   = -( FLX_INT(:,:, 0) + DFDTS (:,:, 0) * DELT)
            if(associated(OLRA )) OLRA  = MAPL_UNDEF
            if(associated(OLC  )) OLC   = -( FLC_INT(:,:, 0) + DFDTSC(:,:, 0) * DELT)
            if(associated(OLA  )) OLA   = MAPL_UNDEF
            if(associated(OLCC5)) then
               where(CLDTT <= 0.05 )
                  OLCC5 = -( FLC_INT(:,:, 0) + DFDTSC(:,:, 0) * DELT)
               elsewhere
                  OLCC5 = MAPL_UNDEF
               endwhere
            endif

            ! fill surface exports

            ! current surface emitted flux derivative wrt surface temperature (+ve)
            ! pmn: should be deprecated ... same as DSFDTS0
            if(associated(DSFDTS)) DSFDTS = -DFDTS(:,:,LM)

            ! surface emitted flux (+ve)
            if(associated(SFCEM)) SFCEM = SFCEM_INT - DFDTS(:,:,LM) * DELT

            ! absorbed (non-reflected) downward surface fluxes
            ! (remember: downward fluxes are not not linearized)
            if(associated(LWS  )) LWS   = FLX_INT(:,:,LM) + SFCEM_INT
            if(associated(LWSA )) LWSA  = MAPL_UNDEF
            if(associated(LCS  )) LCS   = FLC_INT(:,:,LM) + SFCEM_INT
            if(associated(LAS  )) LAS   = MAPL_UNDEF
            if(associated(LCSC5)) then
               where(CLDTT <= 0.05 )
                  LCSC5 = FLC_INT(:,:,LM) + SFCEM_INT
               elsewhere
                  LCSC5 = MAPL_UNDEF
               endwhere
            endif

            ! surface net downward fluxes
            if(associated(FLNS  )) FLNS   = FLX_INT(:,:,LM) + DFDTS (:,:,LM) * DELT
            if(associated(FLNSNA)) FLNSNA = MAPL_UNDEF
            if(associated(FLNSC )) FLNSC  = FLC_INT(:,:,LM) + DFDTSC(:,:,LM) * DELT
            if(associated(FLNSA )) FLNSA  = MAPL_UNDEF

         end if  ! RRTMG

         ! band OLR and/or TBR output
         if ((USE_RRTMG .or. USE_RRTMGP) .and. any_band_output) then

            allocate(OLRB(IM,JM),_STAT)

            if (USE_RRTMGP) then
               call ESMF_UserCompGetInternalState(GC, 'RRTMGP_state', wrap, status)
               VERIFY_(status)
               rrtmgp_state => wrap%ptr
               if (rrtmgp_state%initialized) &
                    band_lims_wvn = rrtmgp_state%k_dist%get_band_lims_wavenumber()
            end if

            do ibnd = 1,nbndlw
               if (band_output(ibnd)) then
                  write(bb,'(I0.2)') ibnd

                  ! get last full calculation
                  call MAPL_GetPointer(INTERNAL, ptr2d, 'OLRB'//bb//'RG', _RC)
                  OLRB = ptr2d

                  ! update for surface temperature on heartbeat
                  call MAPL_GetPointer(INTERNAL, ptr2d, 'DOLRB'//bb//'RGDT', _RC)
                  OLRB = OLRB + ptr2d * DELT

                  ! fill OLRBbbRG if requested
                  call MAPL_GetPointer(EXPORT, ptr2d, 'OLRB'//bb//'RG', _RC)
                  if (associated(ptr2d)) then
                     if (all(OLRB == 0.)) then
                        ! handles pre-first-full-calc case
                        ptr2d = MAPL_UNDEF
                     else
                        ptr2d = OLRB
                     end if
                  end if

                  ! calculate TBRBbbRG if requested
                  call MAPL_GetPointer(EXPORT, ptr2d, 'TBRB'//bb//'RG', _RC)
                  if (associated(ptr2d)) then
                     if (USE_RRTMG) then
                        wn1 = wavenum1(ibnd)*100.; wn2 = wavenum2(ibnd)*100.  ! [m-1]
                        call Tbr_from_band_flux(IM, JM, OLRB, wn1, wn2, ptr2d, _RC)
                     else ! RRTMGP
                        if (rrtmgp_state%initialized) then
                           wn1 = band_lims_wvn(1,ibnd)*100.; wn2 = band_lims_wvn(2,ibnd)*100.  ! [m-1]
                           call Tbr_from_band_flux(IM, JM, OLRB, wn1, wn2, ptr2d, _RC)
                        else
                           ptr2d = MAPL_UNDEF
                        end if
                     end if
                  end if

               end if
            end do

            deallocate(OLRB,_STAT)
         end if

         ! update reference linearization to current temperature
         ! pmn: should be deprecated because its moving along the line passing
         !   through point (TS_INT, SFCEM_INT) with slope -DFDTS (:,:,LM) that
         !   was defined only in the last REFRESH(). Better to just stick with
         !   the exports set in REFRESH() alone.
         if(associated(DSFDTS0)) DSFDTS0 =           - DFDTS(:,:,LM)
         if(associated(SFCEM0 )) SFCEM0  = SFCEM_INT - DFDTS(:,:,LM) * DELT
         if(associated(TSREFF )) TSREFF  = TSINST

         ! Process RAT diagnostics <<>> MSL
         if (nRATS .gt. 0) then
            call MAPL_GetPointer(INTERNAL, DFDTS_RAT,     'DFDTS_RAT',  RC=STATUS)
            call MAPL_GetPointer(INTERNAL, FLX_INT_RAT,   'FLX_RAT',    RC=STATUS)
            call MAPL_GetPointer(INTERNAL, SFCEM_INT_RAT, 'SFCEM_RAT',  RC=STATUS)
            call MAPL_GetPointer(INTERNAL, FLXU_INT_RAT,  'FLXU_RAT',   RC=STATUS)
            call MAPL_GetPointer(IMPORT, EMIS,   'EMIS',   _RC)
            do n=1,nRATS
               ! OLR
               !<<>>         if (MAPL_am_I_root()) then
               !<<>>            write(*,*) '<<>> alloc? ', allocated(nameRATS), ' n: ', n, ' nRATS: ', nRATS
               !<<>>            if (allocated(nameRATS)) write(*,*) '<<>> nameRATS: ', trim(nameRATS(n))
               !<<>>         endif
               gen_str = 'dOLR_'//trim(nameRATS(n)) !nameRATS is the list of active RAT toggles read from AGCM.rc
               call MAPL_GetPointer(EXPORT,   RAT_2d, trim(gen_str),   RC=STATUS) ! Don't verify.
               if (associated(RAT_2d)) then
                  RAT_2d = -( FLX_INT_RAT(:,:,0,n) )
                  RAT_2d = (-( FLX_INT(:,:, 0))) - RAT_2d
                  RAT_2d => null()
               endif
               gen_str = 'dLWS_'//trim(nameRATS(n)) !nameRATS is the list of active RAT toggles read from AGCM.rc
               call MAPL_GetPointer(EXPORT,   RAT_2d, trim(gen_str),   RC=STATUS) ! Don't verify.
               if (associated(RAT_2d)) then
                  RAT_2d   = (FLX_INT(:,:,LM) + SFCEM_INT)-(FLX_INT_RAT(:,:,LM,n) + SFCEM_INT_RAT(:,:,n))
                  RAT_2d => null()
               endif
               gen_str = 'dFLNS_'//trim(nameRATS(n)) !nameRATS is the list of active RAT toggles read from AGCM.rc
               call MAPL_GetPointer(EXPORT,   RAT_2d, trim(gen_str),   RC=STATUS) ! Don't verify.
               if (associated(RAT_2d)) then
                  RAT_2d = (FLX_INT(:,:,LM)) - (FLX_INT_RAT(:,:,LM,n))
                  RAT_2d => null()
               endif
               gen_str = 'dSFCEM_'//trim(nameRATS(n)) !nameRATS is the list of active RAT toggles read from AGCM.rc
               call MAPL_GetPointer(EXPORT,   RAT_2d, trim(gen_str),   RC=STATUS) ! Don't verify.
               if (associated(RAT_2d)) then
                  RAT_2d = SFCEM_INT! - DFDTS(:,:,LM) * DELT
                  RAT_2d = RAT_2d - (SFCEM_INT_RAT(:,:,n))! - DFDTS_RAT(:,:,LM,n) * DELT)
                  RAT_2d => null()
               endif
               gen_str = 'NETTRAP_'//trim(nameRATS(n)) !nameRATS is the list of active RAT toggles read from AGCM.rc
               call MAPL_GetPointer(EXPORT,   RAT_2d, trim(gen_str),   RC=STATUS) ! Don't verify.
               if (associated(RAT_2d)) then
                  RAT_2d = (FLX_INT(:,:,LM)) - &  ! Net DOWNWARD flux
                       (FLX_INT(:,:, 0))
                  RAT_2d = RAT_2d - &
                       ((FLX_INT_RAT(:,:,LM,n)) - & ! Net DOWNWARD flux without RAT at index "n"
                       (FLX_INT_RAT(:,:, 0,n)))
                  RAT_2d => null()
               endif
               gen_str = 'COLTRAP_'//trim(nameRATS(n)) !nameRATS is the list of active RAT toggles read from AGCM.rc
               !         call MAPL_GetPointer(IMPORT, AREA,   'AREA',   RC=STATUS); VERIFY_(STATUS) ! Uncomment for AREA
               call MAPL_GetPointer(EXPORT,   RAT_3d, trim(gen_str),   RC=STATUS) ! Don't verify.
               if (associated(RAT_3d)) then
                  do K = 1, LM
                     RAT_3d(:,:,K) = &
                          (FLX_INT(:,:,K  )) - &
                          (FLX_INT(:,:,K-1))
                     !               RAT_3d(:,:,K) = AREA(:,:) * ( & ! Uncomment this and comment the line below to multiply by area
                     RAT_3d(:,:,K) = ( &    ! Comment this and uncomment the line above to multiply by area
                          RAT_3d(:,:,K) - &
                          ((FLX_INT_RAT(:,:,K  ,n)) - &
                          (FLX_INT_RAT(:,:,K-1,n))))
                  enddo
                  RAT_3d => null()
                  !            AREA => null() ! Uncomment for AREA
               endif
               gen_str = 'FLX_'//trim(nameRATS(n)) !nameRATS is the list of active RAT toggles read from AGCM.rc
               call MAPL_GetPointer(EXPORT,   RAT_3d, trim(gen_str),   RC=STATUS) ! Don't verify.
               if (associated(RAT_3d)) then
                  RAT_3d = FLX_INT(:,:,:) - FLX_INT_RAT(:,:,:,n)
                  RAT_3d => null()
               endif
               gen_str = 'DFDTS_'//trim(nameRATS(n)) !nameRATS is the list of active RAT toggles read from AGCM.rc
               call MAPL_GetPointer(EXPORT,   RAT_3d, trim(gen_str),   RC=STATUS) ! Don't verify.
               if (associated(RAT_3d)) then
                  RAT_3d = DFDTS(:,:,:) - DFDTS_RAT(:,:,:,n)
                  RAT_3d => null()
               endif
               gen_str = 'DELT' !nameRATS is the list of active RAT toggles read from AGCM.rc
               call MAPL_GetPointer(EXPORT,   RAT_2d, trim(gen_str),   RC=STATUS) ! Don't verify.
               if (associated(RAT_2d)) then
                  RAT_2d = DELT
                  RAT_2d => null()
               endif
            enddo
         endif

         !  All done
         deallocate( DUMTT )

         RETURN_(ESMF_SUCCESS)

      end subroutine Update_Flx

   end subroutine RUN

end module GEOS_IrradGridCompMod
