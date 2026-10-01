#include "MAPL.h"

module GEOS_RadiationGridCompMod

   !BOP
   !MODULE: GEOS_RadiationGridCompMod -- Container for atmospheric radiation calculations

   !DESCRIPTION: A composite MAPL/ESMF gridded component (GC) containing the
   !   longwave and shortwave radiation GCs. It is intended as a container
   !   for the ESMF/MAPL Solar and Irrad gridded components used in GEOS.
   !   In SetServices, it creates its children (currently IRRAD and SOLAR;
   !   SATSIM is not yet ported to MAPL3). Its Run method combines results
   !   from the children to produce total radiative exports.
   !
   !   It follows the standard rules for composite ESMF/MAPL GCs. It passes
   !   the ESMF grid that appears in the gridded component to its children,
   !   and all their Imports and Exports are assumed to be on this grid.
   !EOP

   use ESMF
   use MAPL

   use GEOS_IrradGridCompMod, only: irrad_setservices => SetServices
   use GEOS_SolarGridCompMod, only: solarSetServices => SetServices
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
      !   and SOLAR children, connects their flux exports into this
      !   component's own Import state so Run can combine them into the
      !   temperature-tendency exports, and re-exports SOLAR's remaining
      !   exports directly under our own names.
      !EOP

      integer :: status
      integer :: DO_OBIO

      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_INITIALIZE, Initialize, _RC)
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_RUN, Run, _RC)

      ! Create the IRRAD child and invoke its SetServices
      call MAPL_GridCompAddChild(gc, "IRRAD", irrad_setservices, "irrad.yaml", _RC)

      ! Create the SOLAR child and invoke its SetServices
      call MAPL_GridCompAddChild(gc, "SOLAR", solarSetServices, "solar.yaml", _RC)

      ! Pull IRRAD's longwave-flux exports into our own Import state so Run
      ! can combine them into the temperature-tendency exports below.
      call MAPL_GridCompAddConnection(gc, &
           src_comp="IRRAD", &
           src_names="FLX, FLC, FLXA, FLA, DSFDTS0, SFCEM0, TSREFF", &
           dst_comp="<self>", _RC)

      ! Pull SOLAR's shortwave-flux exports into our own Import state so Run
      ! can combine them into the temperature-tendency exports below.
      call MAPL_GridCompAddConnection(gc, &
           src_comp="SOLAR", &
           src_names="FSW, FSC, FSWNA, FSCNA", &
           dst_comp="<self>", _RC)

      ! Re-export SOLAR's other exports directly under our own names, as
      ! the old MAPL2 file's CHILD_ID=SOL-promoted exports used to do.
      call MAPL_GridCompReexport(gc, src_comp="SOLAR", src_name="DRPAR", _RC)
      call MAPL_GridCompReexport(gc, src_comp="SOLAR", src_name="DFPAR", _RC)
      call MAPL_GridCompReexport(gc, src_comp="SOLAR", src_name="DRNIR", _RC)
      call MAPL_GridCompReexport(gc, src_comp="SOLAR", src_name="DFNIR", _RC)
      call MAPL_GridCompReexport(gc, src_comp="SOLAR", src_name="DRUVR", _RC)
      call MAPL_GridCompReexport(gc, src_comp="SOLAR", src_name="DFUVR", _RC)
      call MAPL_GridCompReexport(gc, src_comp="SOLAR", src_name="DRPARN", _RC)
      call MAPL_GridCompReexport(gc, src_comp="SOLAR", src_name="DFPARN", _RC)
      call MAPL_GridCompReexport(gc, src_comp="SOLAR", src_name="DRNIRN", _RC)
      call MAPL_GridCompReexport(gc, src_comp="SOLAR", src_name="DFNIRN", _RC)
      call MAPL_GridCompReexport(gc, src_comp="SOLAR", src_name="DRUVRN", _RC)
      call MAPL_GridCompReexport(gc, src_comp="SOLAR", src_name="DFUVRN", _RC)
      call MAPL_GridCompReexport(gc, src_comp="SOLAR", src_name="FCLD", _RC)
      call MAPL_GridCompReexport(gc, src_comp="SOLAR", src_name="TAUCLI", _RC)
      call MAPL_GridCompReexport(gc, src_comp="SOLAR", src_name="TAUCLW", _RC)
      call MAPL_GridCompReexport(gc, src_comp="SOLAR", src_name="CLDTT", _RC)
      call MAPL_GridCompReexport(gc, src_comp="SOLAR", src_name="ALBEDO", _RC)
      call MAPL_GridCompReexport(gc, src_comp="SOLAR", src_name="FSWBAND", _RC)
      call MAPL_GridCompReexport(gc, src_comp="SOLAR", src_name="FSWBANDNA", _RC)

      ! DROBIO/DFOBIO only exist on SOLAR's side when OBIO coupling is on
      ! (Solar_StateSpecs.rc's COND=SOLAR_TO_OBIO) - mirror the same
      ! resource check here before re-exporting them.
      call MAPL_GridCompGetResource(gc, "USE_OCEANOBIOGEOCHEM", DO_OBIO, default=0, _RC)
      if (DO_OBIO /= 0) then
         call MAPL_GridCompReexport(gc, src_comp="SOLAR", src_name="DROBIO", _RC)
         call MAPL_GridCompReexport(gc, src_comp="SOLAR", src_name="DFOBIO", _RC)
      end if

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
      call MAPL_GridCompGetResource(gc, "ADL_AM1", aam1, default=def_aam1, _RC)
      call MAPL_GridCompGetResource(gc, "ADL_AM2", aam2, default=def_aam2, _RC)
      call MAPL_GridCompGetResource(gc, "ADL_AM30", aam30, default=def_aam30, _RC)
      call MAPL_GridCompGetResource(gc, "ADL_AM4", aam4, default=def_aam4, _RC)
      call MAPL_GridCompGetResource(gc, "RDL_AM1", ram1, default=def_ram1, _RC)
      call MAPL_GridCompGetResource(gc, "RDL_AM2", ram2, default=def_ram2, _RC)
      call MAPL_GridCompGetResource(gc, "RDL_AM30", ram30, default=def_ram30, _RC)
      call MAPL_GridCompGetResource(gc, "RDL_AM4", ram4, default=def_ram4, _RC)
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

      !DESCRIPTION: Runs the IRRAD and SOLAR children and combines their
      !   longwave and shortwave fluxes into pressure-weighted temperature
      !   tendencies.
      !EOP

      integer :: status
      type(ESMF_Grid) :: esmfgrid
      integer :: IM, JM, LM

      ! Pointers to imports (PLEINST plus IRRAD's and SOLAR's exports,
      ! connected into our own Import state in SetServices)
      real, pointer, contiguous, dimension(:, :, :) :: PLE
      real, pointer, contiguous, dimension(:, :, :) :: FLW, FLWCLR, FLWNA, FLA
      real, pointer, contiguous, dimension(:, :, :) :: FSW, FSWCLR, FSWNA, FSCNA
      real, pointer, dimension(:, :) :: DSFDTS, SFCEM, TRD

      ! Scratch pointer for the Edge-array bounds remap below. MAPL always
      ! creates Edge-staggered fields with Fortran bounds 1:LM+1 (no
      ! automatic 0-based remap for VLOC=E fields), but this routine
      ! assumes the 0:LM layer-interface convention throughout (e.g.
      ! PLE(:,:,0:LM-1)/(:,:,1:LM)). A *self*-referencing rank remap
      ! (X(bounds) => X) is rejected by gfortran-15 even when X is
      ! declared CONTIGUOUS ("Rank remapping target must be rank 1 or
      ! simply contiguous") - routing through this distinct
      ! CONTIGUOUS-declared scratch pointer avoids that (verified with an
      ! isolated test compile; see GEOS_IrradGridComp.F90's identical
      ! fix and porting notes for the full explanation).
      real, pointer, contiguous, dimension(:, :, :) :: p3d

      ! Pointers to exports
      real, pointer, dimension(:, :, :) :: RADLW, RADLWC, RADLWNA, RADLWCNA
      real, pointer, dimension(:, :, :) :: RADSW, RADSWC, RADSWNA, RADSWCNA
      real, pointer, dimension(:, :, :) :: DTDT
      real, pointer, dimension(:, :) :: ALW, BLW, RADSRF

      real, allocatable, dimension(:, :, :) :: DMI

      call MAPL_GridCompGet(gc, grid=esmfgrid, num_levels=LM, _RC)
      call MAPL_GridGet(esmfgrid, IM=IM, JM=JM, _RC)

      ! Run the child components
      call MAPL_GridCompTimerStop(gc, "IRRAD", _RC)
      call MAPL_GridCompRunChild(gc, "IRRAD", _RC)
      call MAPL_GridCompTimerStart(gc, "IRRAD", _RC)

      call MAPL_GridCompTimerStop(gc, "SOLAR", _RC)
      call MAPL_GridCompRunChild(gc, "SOLAR", _RC)
      call MAPL_GridCompTimerStart(gc, "SOLAR", _RC)

      ! Get pointers to imports
      call MAPL_StateGetPointer(import, PLE, 'PLEINST', _RC)
      call MAPL_StateGetPointer(import, FLW, 'FLX', _RC)
      call MAPL_StateGetPointer(import, FLWCLR, 'FLC', _RC)
      call MAPL_StateGetPointer(import, FLWNA, 'FLXA', _RC)
      call MAPL_StateGetPointer(import, FLA, 'FLA', _RC)
      call MAPL_StateGetPointer(import, FSW, 'FSW', _RC)
      call MAPL_StateGetPointer(import, FSWCLR, 'FSC', _RC)
      call MAPL_StateGetPointer(import, FSWNA, 'FSWNA', _RC)
      call MAPL_StateGetPointer(import, FSCNA, 'FSCNA', _RC)
      call MAPL_StateGetPointer(import, DSFDTS, 'DSFDTS0', _RC)
      call MAPL_StateGetPointer(import, SFCEM, 'SFCEM0', _RC)
      call MAPL_StateGetPointer(import, TRD, 'TSREFF', _RC)

      ! Edge imports: remap to the 0:LM layer-interface convention used
      ! below (see the declaration comment above for why).
      p3d => PLE;    PLE   (1:IM,1:JM,0:LM) => p3d
      p3d => FLW;    FLW   (1:IM,1:JM,0:LM) => p3d
      p3d => FLWCLR; FLWCLR(1:IM,1:JM,0:LM) => p3d
      p3d => FLWNA;  FLWNA (1:IM,1:JM,0:LM) => p3d
      p3d => FLA;    FLA   (1:IM,1:JM,0:LM) => p3d
      p3d => FSW;    FSW   (1:IM,1:JM,0:LM) => p3d
      p3d => FSWCLR; FSWCLR(1:IM,1:JM,0:LM) => p3d
      p3d => FSWNA;  FSWNA (1:IM,1:JM,0:LM) => p3d
      p3d => FSCNA;  FSCNA (1:IM,1:JM,0:LM) => p3d

      ! Get pointers to exports
      call MAPL_StateGetPointer(export, ALW, 'ALW', _RC)
      call MAPL_StateGetPointer(export, BLW, 'BLW', _RC)
      call MAPL_StateGetPointer(export, RADSRF, 'RADSRF', _RC)
      call MAPL_StateGetPointer(export, DTDT, 'DTDT', _RC)
      call MAPL_StateGetPointer(export, RADLW, 'RADLW', _RC)
      call MAPL_StateGetPointer(export, RADSW, 'RADSW', _RC)
      call MAPL_StateGetPointer(export, RADLWC, 'RADLWC', _RC)
      call MAPL_StateGetPointer(export, RADSWC, 'RADSWC', _RC)
      call MAPL_StateGetPointer(export, RADLWNA, 'RADLWNA', _RC)
      call MAPL_StateGetPointer(export, RADSWNA, 'RADSWNA', _RC)
      call MAPL_StateGetPointer(export, RADLWCNA, 'RADLWCNA', _RC)
      call MAPL_StateGetPointer(export, RADSWCNA, 'RADSWCNA', _RC)

      ! Prepare exports
      if (associated(BLW)) BLW = DSFDTS
      if (associated(ALW)) ALW = SFCEM - DSFDTS * TRD
      if (associated(RADSRF)) RADSRF = FSW(:, :, LM) + FLW(:, :, LM)
      if (associated(DTDT)) DTDT = ((FLW(:, :, 0:LM - 1) - FLW(:, :, 1:LM)) + &
           (FSW(:, :, 0:LM - 1) - FSW(:, :, 1:LM))) * (MAPL_GRAV / MAPL_CP)

      if (associated(RADLW) .or. associated(RADLWC) .or. &
           associated(RADLWNA) .or. associated(RADLWCNA) .or. &
           associated(RADSW) .or. associated(RADSWC) .or. &
           associated(RADSWNA) .or. associated(RADSWCNA)) then

         allocate(DMI(IM, JM, LM), _STAT)
         DMI = MAPL_GRAV / (MAPL_CP * (PLE(:, :, 1:LM) - PLE(:, :, 0:LM - 1)))

         if (associated(RADLW)) RADLW = (FLW(:, :, 0:LM - 1) - FLW(:, :, 1:LM)) * DMI
         if (associated(RADLWC)) RADLWC = (FLWCLR(:, :, 0:LM - 1) - FLWCLR(:, :, 1:LM)) * DMI
         if (associated(RADLWNA)) RADLWNA = (FLWNA(:, :, 0:LM - 1) - FLWNA(:, :, 1:LM)) * DMI
         if (associated(RADLWCNA)) RADLWCNA = (FLA(:, :, 0:LM - 1) - FLA(:, :, 1:LM)) * DMI
         if (associated(RADSW)) RADSW = (FSW(:, :, 0:LM - 1) - FSW(:, :, 1:LM)) * DMI
         if (associated(RADSWC)) RADSWC = (FSWCLR(:, :, 0:LM - 1) - FSWCLR(:, :, 1:LM)) * DMI
         if (associated(RADSWNA)) RADSWNA = (FSWNA(:, :, 0:LM - 1) - FSWNA(:, :, 1:LM)) * DMI
         if (associated(RADSWCNA)) RADSWCNA = (FSCNA(:, :, 0:LM - 1) - FSCNA(:, :, 1:LM)) * DMI

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
