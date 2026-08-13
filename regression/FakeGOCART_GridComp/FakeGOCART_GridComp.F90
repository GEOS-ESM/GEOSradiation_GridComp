#include "MAPL.h"

! Test-double for a real aerosol-optics-provider component (e.g. GOCART2G).
! Publishes a single EXPORT, 'AERO', with itemtype STATE - matching the
! shape of IRRAD's real 'AERO' import - but with the
! 'implements_aerosol_optics_method' Info attribute set to .false., so
! IRRAD's RADIATIVELY_ACTIVE_AEROSOLS branch is skipped without needing a
! real aerosol provider wired into the irrad-sa regression test.

module FakeGOCART_GridCompMod

   use ESMF
   use MAPL, only: MAPL_Verify, MAPL_Return
   use MAPL, only: MAPL_GridCompSetEntryPoint, MAPL_GridCompAddSpec
   use MAPL, only: MAPL_STATEITEM_STATE, MAPL_VERTICAL_STAGGER_CENTER

   implicit none
   private

   public SetServices

contains

   subroutine SetServices(gc, rc)
      type(ESMF_GridComp) :: gc
      integer, intent(out) :: rc

      integer :: status

      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_INITIALIZE, Initialize, _RC)
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_RUN, Run, phase_name="run", _RC)

      call MAPL_GridCompAddSpec(gc, &
           state_intent=ESMF_STATEINTENT_EXPORT, &
           short_name='AERO', &
           standard_name='aerosols', &
           units='kg kg-1', &
           dims='xyz', &
           vertical_stagger=MAPL_VERTICAL_STAGGER_CENTER, &
           itemtype=MAPL_STATEITEM_STATE, _RC)

      _RETURN(_SUCCESS)
   end subroutine SetServices

   subroutine Initialize(gc, import, export, clock, rc)
      type(ESMF_GridComp) :: gc
      type(ESMF_State) :: import
      type(ESMF_State) :: export
      type(ESMF_Clock) :: clock
      integer, intent(out) :: rc

      type(ESMF_State) :: aero
      type(ESMF_Info) :: aero_info
      integer :: status

      call ESMF_StateGet(export, 'AERO', aero, _RC)
      call ESMF_InfoGetFromHost(aero, aero_info, _RC)
      call ESMF_InfoSet(aero_info, key='implements_aerosol_optics_method', value=.false., _RC)

      _RETURN(_SUCCESS)
      _UNUSED_DUMMY(gc)
      _UNUSED_DUMMY(import)
      _UNUSED_DUMMY(clock)
   end subroutine Initialize

   subroutine Run(gc, import, export, clock, rc)
      type(ESMF_GridComp) :: gc
      type(ESMF_State) :: import
      type(ESMF_State) :: export
      type(ESMF_Clock) :: clock
      integer, intent(out) :: rc

      integer :: status

      ! Nothing to do - AERO's Info attribute is set once in Initialize.
      _RETURN(_SUCCESS)
      _UNUSED_DUMMY(gc)
      _UNUSED_DUMMY(import)
      _UNUSED_DUMMY(export)
      _UNUSED_DUMMY(clock)
   end subroutine Run

end module FakeGOCART_GridCompMod

subroutine FakeGOCART_SetServices(gc, rc)
   use ESMF
   use FakeGOCART_GridCompMod, only: mySetServices => SetServices
   type(ESMF_GridComp) :: gc
   integer, intent(out) :: rc
   call mySetServices(gc, rc=rc)
end subroutine FakeGOCART_SetServices
