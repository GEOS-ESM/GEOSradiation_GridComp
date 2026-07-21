#include "MAPL.h"

module GEOS_RadiationGridCompMod

   !BOP
   !MODULE: GEOS_RadiationGridCompMod -- Container for atmospheric radiation calculations

   !DESCRIPTION: A composite MAPL/ESMF gridded component (GC) containing the
   !   longwave and shortwave radiation GCs. It is intended as a container
   !   for the ESMF/MAPL Solar and Irrad gridded components used in GEOS.
   !   In SetServices, it creates its children (currently just IRRAD; SOLAR
   !   and SATSIM are not yet ported to MAPL3). Its Run method combines
   !   results from the children to produce total radiative exports.
   !
   !   It follows the standard rules for composite ESMF/MAPL GCs. It passes
   !   the ESMF grid that appears in the gridded component to its children,
   !   and all their Imports and Exports are assumed to be on this grid.
   !
   !   NOTE: Until SOLAR is ported, only the longwave-only exports (RADLW,
   !   RADLWC, RADLWNA, RADLWCNA, ALW, BLW) are populated by Run. DTDT,
   !   RADSRF, and the shortwave exports (RADSW, RADSWC, RADSWNA, RADSWCNA)
   !   require both longwave and shortwave fluxes and are left unpopulated
   !   until SOLAR is available.
   !EOP

   use ESMF
   use MAPL

   use GEOS_IrradGridCompMod, only: irrad_setservices => SetServices
   ! use GEOS_SolarGridCompMod,  only: solarSetServices  => SetServices  ! NOT ported to MAPL3 yet
   ! use GEOS_SatsimGridCompMod, only: satsimSetServices => SetServices  ! NOT ported to MAPL3 yet

   implicit none
   private

   public SetServices

contains

   !BOP
   !IROUTINE: SetServices -- Sets ESMF services for this component
   !INTERFACE:
   subroutine SetServices(gc, rc)

      !ARGUMENTS:
      type(ESMF_GridComp) :: gc
      integer, intent(out) :: rc

      !DESCRIPTION: Registers Initialize and Run methods, creates the IRRAD
      !   child, and connects IRRAD's longwave-flux exports into this
      !   component's own Import state so Run can combine them into the
      !   temperature-tendency exports.
      !EOP

      integer :: status

      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_INITIALIZE, Initialize, _RC)
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_RUN, Run, _RC)

      ! Create the IRRAD child and invoke its SetServices
      call MAPL_GridCompAddChild(gc, "IRRAD", irrad_setservices, "irrad.yaml", _RC)

      ! Pull IRRAD's longwave-flux exports into our own Import state so Run
      ! can combine them into the temperature-tendency exports below.
      call MAPL_GridCompAddConnection(gc, &
           src_comp="IRRAD", &
           src_names="FLX, FLC, FLXA, FLA, DSFDTS0, SFCEM0, TSREFF", &
           dst_comp="<self>", _RC)

      ! Set the state variable specs generated from Radiation_StateSpecs.rc
#include "Radiation_Import___.h"
#include "Radiation_Export___.h"

      _RETURN(_SUCCESS)

   end subroutine SetServices

   !BOP
   !IROUTINE: Initialize -- Initialize method for the composite Radiation Gridded Component
   !INTERFACE:
   subroutine Initialize(gc, import, export, clock, rc)

      !USES:
      use cloud_condensate_inhomogeneity, only: set_inhomogeneity
      use cloud_subcol_gen, only: initialize_cloud_subcol_gen, &
           def_aam1, def_aam2, def_aam30, def_aam4, &
           def_ram1, def_ram2, def_ram30, def_ram4

      !ARGUMENTS:
      type(ESMF_GridComp) :: gc
      type(ESMF_State) :: import
      type(ESMF_State) :: export
      type(ESMF_Clock) :: clock
      integer, intent(out) :: rc

      !DESCRIPTION: Initializes module-level cloud subcolumn generator
      !   details used by RRTMG[P] LW and SW. This should be done away with
      !   when there is a better treatment of these tables.
      !EOP

      integer :: status

      ! Condensate inhomogeneity type from resource file
      ! ih == 0: homogeneous
      ! ih == 1: inhomogeneous, beta  distribution
      ! ih == 2: inhomogeneous, gamma distribution
      integer :: ih

      ! Correlation length parameters from resource file
      real :: aam1, aam2, aam30, aam4
      real :: ram1, ram2, ram30, ram4

      call MAPL_GridCompGetResource(gc, "RAD_CONDENSATE_INHOMOGENEITY", ih, default=1, _RC)
      call set_inhomogeneity(ih)

      ! Set RRTMG[P] cloud subcolumn generator correlation length parameters
      ! to non-default values from MAPL resource parameters, if given.
      call MAPL_GridCompGetResource(gc, "ADL_AM1",  aam1,  default=def_aam1,  _RC)
      call MAPL_GridCompGetResource(gc, "ADL_AM2",  aam2,  default=def_aam2,  _RC)
      call MAPL_GridCompGetResource(gc, "ADL_AM30", aam30, default=def_aam30, _RC)
      call MAPL_GridCompGetResource(gc, "ADL_AM4",  aam4,  default=def_aam4,  _RC)
      call MAPL_GridCompGetResource(gc, "RDL_AM1",  ram1,  default=def_ram1,  _RC)
      call MAPL_GridCompGetResource(gc, "RDL_AM2",  ram2,  default=def_ram2,  _RC)
      call MAPL_GridCompGetResource(gc, "RDL_AM30", ram30, default=def_ram30, _RC)
      call MAPL_GridCompGetResource(gc, "RDL_AM4",  ram4,  default=def_ram4,  _RC)
      call initialize_cloud_subcol_gen( &
           adl_am1=aam1, adl_am2=aam2, adl_am30=aam30, adl_am4=aam4, &
           rdl_am1=ram1, rdl_am2=ram2, rdl_am30=ram30, rdl_am4=ram4)

      _RETURN(_SUCCESS)
      _UNUSED_DUMMY(import)
      _UNUSED_DUMMY(export)
      _UNUSED_DUMMY(clock)

   end subroutine Initialize

   !BOP
   !IROUTINE: Run -- Run method for the composite Radiation Gridded Component
   !INTERFACE:
   subroutine Run(gc, import, export, clock, rc)

      !ARGUMENTS:
      type(ESMF_GridComp) :: gc
      type(ESMF_State) :: import
      type(ESMF_State) :: export
      type(ESMF_Clock) :: clock
      integer, intent(out) :: rc

      !DESCRIPTION: Runs the IRRAD child and combines its longwave fluxes
      !   into pressure-weighted temperature tendencies. Until SOLAR is
      !   ported, only the longwave-only exports are populated (see module
      !   header).
      !EOP

      integer :: status
      type(ESMF_Grid) :: esmfgrid
      integer :: IM, JM, LM

      ! Pointers to imports (PLEINST plus IRRAD's exports, connected into
      ! our own Import state in SetServices)
      real, pointer, dimension(:,:,:) :: PLE
      real, pointer, dimension(:,:,:) :: FLW, FLWCLR, FLWNA, FLA
      real, pointer, dimension(:,:  ) :: DSFDTS, SFCEM, TRD

      ! Pointers to exports
      real, pointer, dimension(:,:,:) :: RADLW, RADLWC, RADLWNA, RADLWCNA
      real, pointer, dimension(:,:  ) :: ALW, BLW

      real, allocatable, dimension(:,:,:) :: DMI

      call MAPL_GridCompGet(gc, grid=esmfgrid, num_levels=LM, _RC)
      call MAPL_GridGet(esmfgrid, im=IM, jm=JM, _RC)

      ! Run the child components
      call MAPL_GridCompTimerStop(gc, "IRRAD", _RC)
      call MAPL_GridCompRunChild(gc, "IRRAD", _RC)
      call MAPL_GridCompTimerStart(gc, "IRRAD", _RC)

      ! Get pointers to imports
      call MAPL_StateGetPointer(import, PLE, 'PLEINST', _RC)
      call MAPL_StateGetPointer(import, FLW,    'FLX',     _RC)
      call MAPL_StateGetPointer(import, FLWCLR, 'FLC',     _RC)
      call MAPL_StateGetPointer(import, FLWNA,  'FLXA',    _RC)
      call MAPL_StateGetPointer(import, FLA,    'FLA',     _RC)
      call MAPL_StateGetPointer(import, DSFDTS, 'DSFDTS0', _RC)
      call MAPL_StateGetPointer(import, SFCEM,  'SFCEM0',  _RC)
      call MAPL_StateGetPointer(import, TRD,    'TSREFF',  _RC)

      ! Get pointers to exports
      call MAPL_StateGetPointer(export, ALW,      'ALW',      _RC)
      call MAPL_StateGetPointer(export, BLW,      'BLW',      _RC)
      call MAPL_StateGetPointer(export, RADLW,    'RADLW',    _RC)
      call MAPL_StateGetPointer(export, RADLWC,   'RADLWC',   _RC)
      call MAPL_StateGetPointer(export, RADLWNA,  'RADLWNA',  _RC)
      call MAPL_StateGetPointer(export, RADLWCNA, 'RADLWCNA', _RC)

      ! Prepare exports
      if (associated(BLW)) BLW = DSFDTS
      if (associated(ALW)) ALW = SFCEM - DSFDTS*TRD

      if (associated(RADLW) .or. associated(RADLWC) .or. &
          associated(RADLWNA) .or. associated(RADLWCNA)) then

         allocate(DMI(IM,JM,LM), _STAT)
         DMI = MAPL_GRAV/(MAPL_CP*(PLE(:,:,1:LM)-PLE(:,:,0:LM-1)))

         if (associated(RADLW))    RADLW    = (FLW   (:,:,0:LM-1) - FLW   (:,:,1:LM))*DMI
         if (associated(RADLWC))   RADLWC   = (FLWCLR(:,:,0:LM-1) - FLWCLR(:,:,1:LM))*DMI
         if (associated(RADLWNA))  RADLWNA  = (FLWNA (:,:,0:LM-1) - FLWNA (:,:,1:LM))*DMI
         if (associated(RADLWCNA)) RADLWCNA = (FLA   (:,:,0:LM-1) - FLA   (:,:,1:LM))*DMI

         deallocate(DMI, _STAT)

      end if

      call MAPL_GridCompTimerStop(gc, "TOTAL", _RC)

      _RETURN(_SUCCESS)
      _UNUSED_DUMMY(clock)

   end subroutine Run

end module GEOS_RadiationGridCompMod

subroutine Radiation_SetServices(gc, rc)
   use ESMF
   use GEOS_RadiationGridCompMod, only: mySetServices => SetServices
   type(ESMF_GridComp) :: gc
   integer, intent(out) :: rc
   call mySetServices(gc, rc=rc)
end subroutine Radiation_SetServices

