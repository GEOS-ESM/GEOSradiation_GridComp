module rrtmg_lw_rtrnmc

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2009, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! --------- Modules ----------

   use parrrtm, only : nbndlw, ngptlw
   use rrlw_tbl, only: bpade, tblint, tau_tbl, exp_tbl, tfn_tbl
   use rrtmg_lw_setcoef, only : &
     pwvcm, planklay, planklev, plankbnd, dplankbnd_dTs
   use rrlw_con, only: fluxfac
   use rrlw_wvn, only: ngb, delwave

   implicit none 
      
contains

   !-----------------------------------------------------
   subroutine rtrnmc (ncol, nlay, dudTs, &
      semiss, taug, pfracs, cloudy, taucmc, &
      totuflux, totdflux, totuclfl, totdclfl, &
      dtotuflux_dTs, dtotuclfl_dTs, &
      band_output, olrb, dolrb_dTs)
   !-----------------------------------------------------
   !
   !  Original version:   E. J. Mlawer, et al. RRTM_V3.0
   !  Revision for GCMs:  Michael J. Iacono; October, 2002
   !  Revision for F90:  Michael J. Iacono; June, 2006
   !  Revision for dFdT option: M. J. Iacono and E. J. Mlawer, Nov 2009
   !
   !  This program calculates the upward fluxes and downward fluxes
   !  for an arbitrary clear or cloudy atmosphere. The input is the
   !  atmospheric profile, all Planck function information, and the
   !  McICA cloud mask. A variable diffusivity angle (secdiff)
   !  is used for the angle integration. Bands 2-3 and 5-9 use a
   !  value for SECDIFF that varies from 1.50 to 1.80 as a function 
   !  of the column water vapor, and other bands use a value of 1.66.
   !  The Gaussian weight appropriate to this angle (wtdiff=0.5) is
   !  applied here. Note that use of the emissivity angle for the
   !  flux integration can cause errors of 1 to 4 W/m2 within cloudy
   !  layers. Clouds are treated with the McICA stochastic approach.
   !
   !  Also provides the optional capability to calculate the deriv-
   !  ative of upward flux respect to surface temperature using the
   !  pre-tabulated derivative of the Planck function with respect 
   !  to temperature integrated over each spectral band.
   ! 
   !  PMN 2021/06 This version depends explicity on binary subcolumn
   !  clouds i.e., cldfmc in {0,1}. Thus, no clear / cloudy portion
   !  breakdown is needed and the optical cloudiness of a subcolumn
   !  / g-point gridcell is totally specified by taucmc = 0. or > 0.
   !-------------------------------------------------------

      integer, intent(in) :: ncol   ! number of columns
      integer, intent(in) :: nlay   ! number of layers
      logical, intent(in) :: dudTs  ! do d(upflux)/d(Tsurf) calcs

      real,    intent(in) :: semiss  (     nbndlw,ncol)  ! surface emissivity
      real,    intent(in) :: taug    (nlay,ngptlw,ncol)  ! gas optical depth
      real,    intent(in) :: pfracs  (nlay,ngptlw,ncol)  ! Planck fractions
      logical, intent(in) :: cloudy  (nlay,       ncol)  ! cloudy for ANY g-point
      real,    intent(in) :: taucmc  (nlay,ngptlw,ncol)  ! cloud optical thickness
     
      ! spectrally summed fluxes and upward flux derivatives wrt Tsurf
      real, intent(out) :: totuflux      (0:nlay,ncol)  ! upward longwave flux [W/m2]
      real, intent(out) :: totdflux      (0:nlay,ncol)  ! downward longwave flux [W/m2]
      real, intent(out) :: totuclfl      (0:nlay,ncol)  ! clrsky upward lw flux [W/m2]
      real, intent(out) :: totdclfl      (0:nlay,ncol)  ! clrsky downward lw flux [W/m2]
      real, intent(out) :: dtotuflux_dTs (0:nlay,ncol)  ! d/d(Tsurf) [W/m2/K]
      real, intent(out) :: dtotuclfl_dTs (0:nlay,ncol)  ! d/d(Tsurf) [W/m2/K]

      ! which band OLRs to calculate?
      logical, intent(in) :: band_output (nbndlw)

      ! band OLRs and d/dTs
      real, intent(out) :: olrb      (nbndlw,ncol)  ! [W/m2]
      real, intent(out) :: dolrb_dTs (nbndlw,ncol)  ! [W/m2/K]

      ! ----- Local -----
   
      real :: agas(nlay)    ! gas absorptivity
      real :: atot(nlay)    ! mixed gas plus cloud absorptivity
      real :: bbugas(nlay)  ! gas Planck function for upward rt
      real :: bbutot(nlay)  ! gas and cloud Planck function for upward rt
     
      real :: secdiff    ! diffusivity angle
      real :: odepth     ! gas optical depth
      real :: odcld      ! cloud optical depth
      real :: odtot      ! optical depth of gas and cloud
      real :: tfacgas    ! gas pade factor, used for planck fn
      real :: tfactot    ! gas and cloud pade factor, used for planck fn
      real :: bbdgas     ! gas planck function for downward rt
      real :: gassrc     ! source radiance due to gas only
      real :: bbdtot     ! gas and cloud planck function for downward rt
      real :: tblind     ! real lookup table index
      real :: radld      ! downward radiance
      real :: radclrd    ! downward radiance for clear column
      real :: radlu      ! upward radiance
      real :: radclru    ! upward radiance for clear column
      real :: rad0       ! surface emitted radiance
      real :: reflect    ! surface reflectance

      real :: plfrac, blay, dplankup, dplankdn
      real :: sumfac, deluflux, deluderiv

      ! derivatives with respect to surface temperature
      real :: d_rad0_dTs, d_radlu_dTs, d_radclru_dTs

      integer :: icol          ! column index
      integer :: lev           ! level index
      integer :: ig            ! g-point index
      integer :: ibnd          ! band index
      integer :: ittot, itgas  ! lookup table indices
   
      ! flags when clear-/total-sky downward streams diverge
      logical :: down_streams_diverge

      ! This weight corresponds to the standard sec(diffusivity angle) of 1.66.
      real, parameter :: wtdiff = 0.5  ! weight for radiance to flux conversion

      ! diffusivity factors
      real, parameter, dimension(16) :: a0 = &
         [1.66 , 1.55 , 1.58 , 1.66 , &
          1.54 , 1.454, 1.89 , 1.33 , &
          1.668, 1.66 , 1.66 , 1.66 , &
          1.66 , 1.66 , 1.66 , 1.66 ]
      real, parameter, dimension(16) :: a1 = &
         [0.00 , 0.25 , 0.22 , 0.00 , &
          0.13 , 0.446,-0.10 , 0.40 , &
         -0.006, 0.00 , 0.00 , 0.00 , &
          0.00 , 0.00 , 0.00 , 0.00 ]
      real, parameter, dimension(16) :: a2 = &
         [0.00 ,-12.0 ,-11.7 , 0.00 , &
         -0.72 ,-0.243, 0.19 ,-0.062, &
          0.414, 0.00 , 0.00 , 0.00 , &
          0.00 , 0.00 , 0.00 , 0.00 ]

      ! zero spectral sum accumulators
      totuflux = 0.
      totdflux = 0.
      totuclfl = 0.
      totdclfl = 0.
      if (dudTs) then
         dtotuflux_dTs = 0.
         dtotuclfl_dTs = 0.
      end if

      ! zero band output accumulators
      if (any(band_output)) then
         olrb = 0.
         if (dudTs) dolrb_dTs = 0.
      end if 

      ! main loop
      ! g-points are "monochromatic" in CKD
      do icol = 1,ncol
      do ig = 1,ngptlw 

         ibnd = ngb(ig)

         ! factor for accumulating fluxes across g-points
         sumfac = wtdiff * delwave(ibnd) * fluxfac
      
         ! Reset diffusivity angle for Bands 2-3 and 5-9 to vary (between 1.50
         ! and 1.80) as a function of total column water vapor. The function
         ! has been defined to minimize flux and cooling rate errors in these
         ! bands over a wide range of precipitable water values.

         if (ibnd == 1 .or. ibnd == 4 .or. ibnd >= 10) then
            secdiff = 1.66 
         else
            secdiff = a0(ibnd) + a1(ibnd) * exp(a2(ibnd)*pwvcm(icol))
            if (secdiff > 1.80) then
               secdiff = 1.80 
            else if (secdiff < 1.50) then
               secdiff = 1.50 
            end if
         endif
    
         ! TOA downward LW is zero
         ! and column starts "clear"
         radld = 0.; radclrd = 0. 

         ! see below ... this flag is set false at the top of every SUB-column
         ! (g-point) of the current gridcolumn, but becomes true at the first
         ! layer (going down) in which ANY of the sub-columns is cloudy.
         down_streams_diverge = .false.

         ! Downward radiative transfer loop
         do lev = nlay,1,-1

            plfrac = pfracs(lev,ig,icol)
            blay = planklay(ibnd,lev,icol)
            dplankup = planklev(ibnd,lev,  icol) - blay
            dplankdn = planklev(ibnd,lev-1,icol) - blay

            odepth = secdiff * taug(lev,ig,icol)
            if (odepth < 0.) odepth = 0. 

            ! ---------------------------------------------------
            ! Table lookup for efficiency based on optical depth:
            ! ---------------------------------------------------
            ! Optical depth, tau, is converted to a discretized
            ! table index, itau in [0:ntbl], by
            !    tblind = tau / (bpade + tau)
            !    itau = tblint * tblind + 0.5
            ! where bpade = 1/(pade constant),
            !       tblint = float(ntbl),
            !   and [0:ntbl] is the dimension of the tables.
            !
            ! The three look-up-tables (LUTs) are:
            ! (1) tau_tbl = Optical depth LUT.
            !        Internally
            !           tau_tbl(itau) = bpade * chi/(1-chi)
            !        where chi = itau / tblint.
            !     As such, tau_tbl(itau) is a discretized
            !     version of tau consistent with itau. Since
            !     it is this descretized tau that is used for
            !     the other two table generations (see below)
            !     it is important to use this descretized tau
            !     externally (rather than the original tau),
            !     for consistency.
            !
            ! (2) exp_tbl = Exponential LUT for transmittance.
            !        Internally = exp(-tau_tbl(itau)).
            !
            ! (3) tfn_tbl = Tau transition function LUT.
            !     (i.e. the transition of the Planck function
            !     from that for the mean layer temperature to
            !     that for the layer boundary temperature as
            !     a function of optical depth. The "linear in
            !     tau" method is used to make the table.)
            !     Internally also a function of tau_tbl(itau).
            ! ---------------------------------------------------

            ! gas only
            tblind = odepth / (bpade + odepth)
            itgas = tblint * tblind + 0.5 
            agas(lev) = 1. - exp_tbl(itgas)
            tfacgas = tfn_tbl(itgas)
            bbdgas      = plfrac * (blay + tfacgas * dplankdn)
            bbugas(lev) = plfrac * (blay + tfacgas * dplankup)

            ! Update of LW, rad, through a layer of transmissivity T:
            !   rad -> rad * T + bbd * (1-T) = rad + (bbd - rad) * A
            ! where bbd is effective Planck func and A = 1-T.

            if (taucmc(lev,ig,icol) <= 0.) then

               ! through a clear layer, so gas only
               radld = radld + (bbdgas - radld) * agas(lev)

            else

               ! through a cloudy layer, so mixed gas plus cloud ...
               ! REMEMBER: binary clouds, no clear/cloudy portion

               odcld = secdiff * taucmc(lev,ig,icol)

               ! since cloud optical depth must be added in, we should
               ! add it to the *discretized* gas optical depth, since
               ! the latter is what is actually used to generate the
               ! exp_tbl and tfn_tbl table values above.
               odepth = tau_tbl(itgas)
               odtot = odepth + odcld

               tblind = odtot / (bpade + odtot)
               ittot = tblint * tblind + 0.5 
               atot(lev) = 1. - exp_tbl(ittot)
               tfactot = tfn_tbl(ittot)
               bbdtot      = plfrac * (blay + tfactot * dplankdn)
               bbutot(lev) = plfrac * (blay + tfactot * dplankup)

               radld = radld + (bbdtot - radld) * atot(lev)

            endif

            ! spectrally integrate downward flux (total-sky)
            totdflux(lev-1,icol) = totdflux(lev-1,icol) + sumfac * radld

            ! total-sky and clear-sky fluxes (for downward, totdflux and
            ! totdclfl, respectively) are spectrally integrated across
            ! all gpoints (subcolumns). Therefore, the downward total-
            ! and clear-sky flux streams will be the same at the TOA and
            ! working down through completely clear layers. But once a
            ! layer is reached which is cloudy in at least one subcolumn,
            ! the two streams diverge for that layer and all layers below.

            if (.not. down_streams_diverge) then
               if (cloudy(lev,icol)) down_streams_diverge = .true.
            end if
            if (down_streams_diverge) then
               radclrd = radclrd + (bbdgas - radclrd) * agas(lev) 
            else
               radclrd = radld
            endif

            ! spectrally integrate downward flux (clear-sky)
            totdclfl(lev-1,icol) = totdclfl(lev-1,icol) + sumfac * radclrd

         enddo  ! downward loop

         ! Include the contribution of spectrally varying longwave emissivity
         ! and reflection from the surface to the upward radiative transfer.
         ! Note: Spectral and Lambertian reflection are identical for the
         ! diffusivity angle flux integration used here.
         ! Note: The emissivity is applied to plankbnd and dplankbnd_dTs when 
         ! they are defined in subroutine setcoef. 
    
         ! surface emission
         rad0 = pfracs(1,ig,icol) * plankbnd(ibnd,icol)
         if (dudTs) d_rad0_dTs = pfracs(1,ig,icol) * dplankbnd_dTs(ibnd,icol)

         ! Add in specular reflection of surface downward radiance
         reflect = 1. - semiss(ibnd,icol)
         radlu   = rad0 + reflect * radld
         radclru = rad0 + reflect * radclrd
         totuflux(0,icol) = totuflux(0,icol) + sumfac * radlu
         totuclfl(0,icol) = totuclfl(0,icol) + sumfac * radclru
         if (dudTs) then
            d_radlu_dTs   = d_rad0_dTs
            d_radclru_dTs = d_rad0_dTs
            dtotuflux_dTs(0,icol) = dtotuflux_dTs(0,icol) + sumfac * d_radlu_dTs
            dtotuclfl_dTs(0,icol) = dtotuclfl_dTs(0,icol) + sumfac * d_radclru_dTs
         endif

         ! Upward radiative transfer loop
         do lev = 1,nlay

            if (taucmc(lev,ig,icol) <= 0.) then

               ! clear subcolumn layer
               radlu = radlu + (bbugas(lev) - radlu) * agas(lev)
               if (dudTs) d_radlu_dTs = d_radlu_dTs - d_radlu_dTs * agas(lev)

            else

               ! cloudy subcolumn layer
               radlu = radlu + (bbutot(lev) - radlu) * atot(lev)
               if (dudTs) d_radlu_dTs = d_radlu_dTs - d_radlu_dTs * atot(lev)

            endif

            ! spectrally integrate upward flux (total-sky)
            deluflux = sumfac * radlu
            totuflux(lev,icol) = totuflux(lev,icol) + deluflux

            ! If downward streams diverged, then upward streams must also
            ! because surface reflected flux is different for each stream.

            if (down_streams_diverge) then
               radclru = radclru + (bbugas(lev) - radclru) * agas(lev) 
            else
               radclru = radlu
            endif

            ! spectrally integrate upward flux (clear-sky)
            totuclfl(lev,icol) = totuclfl(lev,icol) + sumfac * radclru

            if (dudTs) then
               if (down_streams_diverge) then
                  d_radclru_dTs = d_radclru_dTs - d_radclru_dTs * agas(lev)
               else
                  d_radclru_dTs = d_radlu_dTs
               endif
               deluderiv = sumfac * d_radlu_dTs
               dtotuflux_dTs(lev,icol) = dtotuflux_dTs(lev,icol) + deluderiv
               dtotuclfl_dTs(lev,icol) = dtotuclfl_dTs(lev,icol) + sumfac * d_radclru_dTs
            endif

         enddo  ! layer
          
         ! OLR band output (we are at TOA)
         if (band_output(ibnd)) then
            olrb(ibnd,icol) = olrb(ibnd,icol) + deluflux
            if (dudTs) dolrb_dTs(ibnd,icol) = dolrb_dTs(ibnd,icol) + deluderiv
         end if

      end do  ! g-point
      end do  ! column

   end subroutine rtrnmc

end module rrtmg_lw_rtrnmc
