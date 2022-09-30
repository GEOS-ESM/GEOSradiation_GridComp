!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2009, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------
!
! ****************************************************************************
! *                                                                          *
! *                             RRTMG_SW                                     *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                 a rapid radiative transfer model                         *
! *                  for the solar spectral region                           *
! *           for application to general circulation models                  *
! *                                                                          *
! *                                                                          *
! *           Atmospheric and Environmental Research, Inc.                   *
! *                       131 Hartwell Avenue                                *
! *                       Lexington, MA 02421                                *
! *                                                                          *
! *                                                                          *
! *                          Eli J. Mlawer                                   *
! *                       Jennifer S. Delamere                               *
! *                        Michael J. Iacono                                 *
! *                        Shepard A. Clough                                 *
! *                       David M. Berthiaume                                *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                      email:  miacono@aer.com                             *
! *                      email:  emlawer@aer.com                             *
! *                      email:  jdelamer@aer.com                            *
! *                                                                          *
! *       The authors wish to acknowledge the contributions of the           *
! *       following people:  Steven J. Taubman, Patrick D. Brown,            *
! *       Ronald E. Farren, Luke Chen, Robert Bergstrom.                     *
! *                                                                          *
! ****************************************************************************

#include "MAPL_Generic.h"
    
module rrtmg_sw_rad

   use ESMF
   use MAPL

   use rrsw_vsn
   use cloud_subcol_gen, only: &
      generate_stochastic_clouds, clearCounts_threeBand
   use rrtmg_sw_cldprmc, only: cldprmc_sw
   use rrtmg_sw_spcvmc, only: spcvmc_sw
   use iso_fortran_env, only: error_unit

   implicit none

   public :: rrtmg_sw, earth_sun

contains

   subroutine rrtmg_sw (MAPL, &
      rpart, ncol, nlay, &
      scon, adjes, coszen, isolvar, &
      play, plev, tlay, &
      h2ovmr, o3vmr, co2vmr, ch4vmr, o2vmr, &
      iceflgsw, liqflgsw, &
      cld, ciwp, clwp, rei, rel, &
      dyofyr, zm, alat, &
      iaer, tauaer, ssaaer, asmaer, &
      asdir, asdif, aldir, aldif, &
      cloudLM, cloudMH, normFlx, &
      clearCounts, swuflx, swdflx, swuflxc, swdflxc, &
      nirr, nirf, parr, parf, uvrr, uvrf, &
      tautp, tauhp, taump, taulp, &
      do_FAR, taumol_age, taumol_age_limit, &
      taur, taug, sflxzen, ssi, &
      taucld_age, taucld_age_limit, &
      Rcldycol, Rcldymc, taucmc, ssacmc, asmcmc, taormc, &
      bndscl, indsolvar, solcycfrac, &  ! optional inputs
      RC)

      use parrrsw, only : nbndsw, ngptsw

      ! ----- Inputs -----

      type(MAPL_MetaComp), intent(inout) :: MAPL

      ! dimensions
      ! ----------
      integer, intent(in) :: rpart                   ! Number of columns in a partition
      integer, intent(in) :: ncol                    ! Number of horizontal columns     
      integer, intent(in) :: nlay                    ! Number of model layers

      ! orbit
      ! -----
      real, intent(in) :: scon                       ! Solar constant (W/m2)

      ! Notes: SCON is total solar irradiance averaged over the solar cycle. If scon == 0,
      ! an internal solar constant will be used. This depends on the value of isolvar. For
      ! isolvar=-1, scon=1368.22 Wm-2 (Kurucz), while for isolvar>=0, scon=1360.85 Wm-2
      ! (NRLSSI2). If scon > 0.0, this internal solar constant will be scaled to the value
      ! of SCON provided.

      real, intent(in) :: adjes                      ! Flux adjustment for Earth/Sun distance
                                                     !    i.e., ~ 1/dist(Earth-Sun in AU)^2
      real, intent(in) :: coszen (ncol)              ! Cosine of solar zenith angle

      ! solar variability
      ! -----------------
      integer, intent(in) :: isolvar                 ! Flag for solar variability method

      ! Notes: isolvar = -1 uses the Kurucz source function, while isolvar >= 0 use the NRLSSI2 
      ! solar model. First, the behavior for SCON = 0: isolvar in {-1,0,3} all have a fixed solar
      ! input at 1AU. For isolvar = -1 it is the Kurucz solar constant of 1368.22 Wm-2, while for
      ! isolvar = {0,3} it is the NRLSSI2 solar constant of 1360.85 Wm-2 (for the 100-50000 cm-1
      ! spectral range only, and based on the mean of solar cycle 13-24). The spectral breakdown
      ! is hardcoded into the model, including the "quiet sun, faculae brightening and sunspot
      ! darkening terms" for NRLSSI2, which are based on the means over solar cycles 13-24. But
      ! an optional scaling by band (via multiplier bndscl(:)) can be applied for isolvar {-1,3}.
      ! For SCON > 0, the behavior of these three cases (isolvar {-1,0,3}) is simple: the solar
      ! input is just scaled to the provided SCON, uniformly across the spectrum (and for each of
      ! the quiet sun, faculae and sunpot terms for NRLSSI2). This scaling is applied before the
      ! optional bndscl scaling.
      !    The isolvar = {1,2} cases require more detailed explanation. Both use the NRLSSI2 model,
      ! and but allow the mix of the quiet sun, faculae, and sunspot terms to vary with time. Since
      ! each of these terms has a different spectral response, this gives a solar input spectra that
      ! varies with time. This is achieved by setting the factors svar_{i,f,s}, which are multipliers
      ! to the average spectral response of the quiet sun, faculae, and sunspot terms (with integral
      ! values {I,F,S}int W/m2.)
      !    The isolvar = 1 case is designed to allow simulations over a generic solar cycle. This
      ! is an 11-year long average solar cycle derived from actual solar cycles 13-24. We denote
      ! this cycle "AvgCyc11". This cycle is implemented as follows: an 11-year cycle is provided
      ! for the faculae and sunpot indices Mg and SB via the {mg,sb}avgcyc arrays in the NRLSSI2
      ! module, which tabulate the index values per month for the 11 years. In conjunction with
      ! this, a linear relationship is provided between Mg and svar_f and between SB and svar_s.
      ! The main input for isolvar = 1 is solcycfrac in [0,1], which is the normalized input
      ! position in the 11-year average cycle. This solcycfrac is used to interpolate into the
      ! {mg,sb}avgcyc arrays to values Mg and SB, which are then converted to svar_f and svar_s
      ! via the linear relationship provided. These linear relationships are designed such that
      ! svar_f = 1 at Mg = <Mg>, the time average Mg index over AvgCyc11, and similarly for svar_s
      ! and SB, such that <svar_{s,f}> are both unity. So, if solcyclfr is uniformly cycled in
      ! [0,1] by the caller of rrtmg_sw_rad(), then each of the faculae and sunspot terms will
      ! also cycle, providing a time and spectrally varying solar cycle, but still with average
      ! flux contributions over the cycle of Fint and Sint (for SCON.eq.0). An extra facility is
      ! provided via the optional indsolar(2) array, which allows for a further multiplication
      ! of the svar_{f,s} terms after they are formed by the linearization above. This multiplier
      ! is time varying, being one at solar minimum, and the value of indsolvar(1,2) at solar
      ! maximum (for 1=Mg, 2=SB respecively), and to vary linearly with solcycfrac between those
      ! extrema (see NRLSSI2 module for further details).
      !    Still discussing isolvar = 1, for SCON.eq.0 we take the hint from not explicitly
      ! setting SCON to let an indsolvar.ne.1 choice cause a deviation from the internal solar 
      ! constant, because, while the svar_{f,s} average to unity over a cycle *without* a time-
      ! varying indsolvar multiplier, they do not do so with it. If, on the other hand, a value
      ! SCON > 0 is provided, we ASSUME that it is a value that we should honor as a MEAN over
      ! the AvgCyc11 cycle. We do this by adjusting a time-invariant quiet sun svar_i factor
      ! to honor the provided mean SCON, even in the presence of indsolvar.ne.1. Further
      ! details are found in the code.
      !    Finally, the insolvar = 2 case. This is a data-driven case provided for when values
      ! of TSI (SCON) and Mg and SB indices are available from data. In this case, the AvgCyc11
      ! cycle is not used explicitly, but the index-to-svar linearizations ARE used. The input
      ! indsolvar in this case provides the (1=Mg,2=SB) indices and produces concomitant mult-
      ! ipliers svar_{f,s} by these linear relationships. Again, for SCON.eq.0 we take the hint
      ! from not explicitly setting SCON to let the solar constant vary naturally as a result
      ! of the above data driven svar_{f,s} values (and by keeping svar_i = 1). BUT, for input
      ! SCON > 0, we take this as the TIME-VARYING, DATA-SUPPLIED TSI and honor it at every time
      ! (not just in a mean sense as for isolvar = 1). We do this by adjusting the quiet sun
      ! svar_i term (see the code).

      real, intent(in), optional :: indsolvar (2)    ! Facular and sunspot amplitude scale facs (isolvar=1),
                                                     !    or Mg and SB indices (isolvar=2)
      real, intent(in), optional :: bndscl (nbndsw)  ! Scale factors for each band
      real, intent(in), optional :: solcycfrac       ! Fraction of averaged 11-year solar cycle (0-1)
                                                     !    at current time (isolvar=1)
                                                     !    0. represents the first day of year 1
                                                     !    1. represents the last day of year 11

      ! profile
      ! -------
      real, intent(in) :: play   (ncol,nlay)         ! Layer pressures (hPa)
      real, intent(in) :: plev   (ncol,nlay+1)       ! Interface pressures (hPa)
      real, intent(in) :: tlay   (ncol,nlay)         ! Layer temperatures (K)

      ! gases
      ! -----
      real, intent(in) :: h2ovmr (ncol,nlay)         ! H2O volume mixing ratio
      real, intent(in) :: o3vmr  (ncol,nlay)         ! O3 volume mixing ratio
      real, intent(in) :: co2vmr (ncol,nlay)         ! CO2 volume mixing ratio
      real, intent(in) :: ch4vmr (ncol,nlay)         ! Methane volume mixing ratio
      real, intent(in) :: o2vmr  (ncol,nlay)         ! Oxygen volume mixing ratio

      ! cloud optics flags
      ! ------------------
      integer, intent(in) :: iceflgsw                ! Flag for ice particle specification
      integer, intent(in) :: liqflgsw                ! Flag for liquid droplet specification

      ! clouds
      ! ------
      real, intent(in) :: cld    (ncol,nlay)         ! Cloud fraction
      real, intent(in) :: ciwp   (ncol,nlay)         ! In-cloud ice water path (g/m2)
      real, intent(in) :: clwp   (ncol,nlay)         ! In-cloud liquid water path (g/m2)
      real, intent(in) :: rei    (ncol,nlay)         ! Cloud ice effective radius (microns)
      real, intent(in) :: rel    (ncol,nlay)         ! Cloud water drop effective radius (microns)

      ! cloud overlap (exponential)
      ! ---------------------------
      integer, intent(in) :: dyofyr                  ! Day of the year
      real, intent(in) :: zm     (ncol,nlay)         ! Heights of level midpoints
      real, intent(in) :: alat   (ncol)              ! Latitude of column [radians]

      ! aerosols (optical props, non-delta-scaled)
      ! ------------------------------------------
      integer, intent(in) :: iaer                    ! aerosol flag (0=off, 10=on)
      real, intent(in) :: tauaer (ncol,nlay,nbndsw)  ! aer optical depth      (iaer=10 only)
      real, intent(in) :: ssaaer (ncol,nlay,nbndsw)  ! aer single scat albedo (iaer=10 only)
      real, intent(in) :: asmaer (ncol,nlay,nbndsw)  ! aer asymmetry param    (iaer=10 only)

      ! surface albedos
      ! ---------------
      real, intent(in) :: asdir  (ncol)              ! UV/vis  surface albedo: direct rad
      real, intent(in) :: asdif  (ncol)              ! UV/vis  surface albedo: diffuse rad
      real, intent(in) :: aldir  (ncol)              ! Near-IR surface albedo: direct rad
      real, intent(in) :: aldif  (ncol)              ! Near-IR surface albedo: diffuse rad

      ! etc
      ! ---
      ! pressure super-layer interface levels for cloud fractions
      integer, intent(in) :: cloudLM  ! Low-mid
      integer, intent(in) :: cloudMH  ! Mid-high

      integer, intent(in) :: normFlx                 ! Normalize fluxes?
                                                     !   0 = no normalization
                                                     !   1 = normalize (by scon*coszen)
      ! FAR controls
      logical, intent(in) :: do_FAR
      real, intent(in) :: taumol_age_limit
      real, intent(in) :: taucld_age_limit

      ! ----- Outputs -----

      ! subcolumn clear counts for Tot|High|Mid|Low super-layers
      integer, intent(out) :: clearCounts(ncol,4)

      real, intent(out) :: swuflx  (ncol,nlay+1)     !   All-sky SW up   flux (W/m2)
      real, intent(out) :: swdflx  (ncol,nlay+1)     !   All-sky SW down flux (W/m2)
      real, intent(out) :: swuflxc (ncol,nlay+1)     ! Clear-sky SW up   flux (W/m2)
      real, intent(out) :: swdflxc (ncol,nlay+1)     ! Clear-sky SW down flux (W/m2)

      ! Output added for Land/Surface process (all-sky)
      real, intent(out) :: nirr    (ncol)            ! Near-IR direct  down SW flux (W/m2)
      real, intent(out) :: nirf    (ncol)            ! Near-IR diffuse down SW flux (W/m2)
      real, intent(out) :: parr    (ncol)            ! Visible direct  down SW flux (W/m2)
      real, intent(out) :: parf    (ncol)            ! Visible diffuse down SW flux (W/m2)
      real, intent(out) :: uvrr    (ncol)            ! UV      direct  down SW flux (W/m2)
      real, intent(out) :: uvrf    (ncol)            ! UV      diffuse down SW flux (W/m2)

      ! In-cloud PAR optical thickness for Tot|High|Mid|Low super-layers
      real, intent(out), dimension (ncol) :: tautp, tauhp, taump, taulp

      integer, intent(out), optional :: RC  ! return code

      ! ------- FAR InOuts -------
      ! if (.not.do_FAR) these can be unassociated pointers since not used

      real, intent(inout), dimension(:),     pointer :: taumol_age      !             (ncol) if (do_FAR)
      real, intent(inout), dimension(:,:,:), pointer :: taur, taug      ! (nlay,ngptsw,ncol) if (do_FAR)
      real, intent(inout), dimension(:,:),   pointer :: sflxzen, ssi    ! (     ngptsw,ncol) if (do_FAR)
      real, intent(inout), dimension(:),     pointer :: taucld_age      !             (ncol) if (do_FAR)
      REAL, intent(inout), dimension(:),     pointer :: Rcldycol        !             (ncol) if (do_FAR)
      REAL, intent(inout), dimension(:,:,:), pointer :: Rcldymc         ! (nlay,ngptsw,ncol) if (do_FAR)
      real, intent(inout), dimension(:,:,:), pointer :: taucmc, ssacmc  ! (nlay,ngptsw,ncol) if (do_FAR)
      real, intent(inout), dimension(:,:,:), pointer :: asmcmc, taormc  ! (nlay,ngptsw,ncol) if (do_FAR)

      ! ----- Locals -----

      integer :: pncol
      integer :: STATUS  ! for MAPL error reporting
      
      ! ASSERTs to catch unphysical or invalid inputs
      _ASSERT(all(play   >= 0.), 'negative values in input:   play')
      _ASSERT(all(plev   >= 0.), 'negative values in input:   plev')
      _ASSERT(all(tlay   >= 0.), 'negative values in input:   tlay')
      _ASSERT(all(h2ovmr >= 0.), 'negative values in input: h2ovmr')
      _ASSERT(all( o3vmr >= 0.), 'negative values in input:  o3vmr')
      _ASSERT(all(co2vmr >= 0.), 'negative values in input: co2vmr')
      _ASSERT(all(ch4vmr >= 0.), 'negative values in input: ch4vmr')
      _ASSERT(all( o2vmr >= 0.), 'negative values in input:  o2vmr')
      _ASSERT(all( asdir >= 0.), 'negative values in input:  asdir')
      _ASSERT(all( aldir >= 0.), 'negative values in input:  aldir')
      _ASSERT(all( asdif >= 0.), 'negative values in input:  asdif')
      _ASSERT(all( aldif >= 0.), 'negative values in input:  aldif')
      _ASSERT(all(   cld >= 0.), 'negative values in input:    cld')
      _ASSERT(all(  ciwp >= 0.), 'negative values in input:   ciwp')
      _ASSERT(all(  clwp >= 0.), 'negative values in input:   clwp')
      _ASSERT(all(   rei >= 0.), 'negative values in input:    rei')
      _ASSERT(all(   rel >= 0.), 'negative values in input:    rel')
      _ASSERT(all(tauaer >= 0.), 'negative values in input: tauaer')
      _ASSERT(all(ssaaer >= 0.), 'negative values in input: ssaaer')

      ! check FAR inputs
      if (do_FAR) then
         _ASSERT(associated(taumol_age),'not associated when do_FAR: taumol_age')
         _ASSERT(all(shape(taumol_age) == [ncol]),'mal-dimensioned: taumol_age')
         _ASSERT(associated(taur),'not associated when do_FAR: taur')
         _ASSERT(all(shape(taur) == [nlay,ngptsw,ncol]),'mal-dimensioned: taur')
         _ASSERT(associated(taug),'not associated when do_FAR: taug')
         _ASSERT(all(shape(taug) == [nlay,ngptsw,ncol]),'mal-dimensioned: taug')
         _ASSERT(associated(sflxzen),'not associated when do_FAR: sflxzen')
         _ASSERT(all(shape(sflxzen) == [ngptsw,ncol]),'mal-dimensioned: sflxzen')
         _ASSERT(associated(ssi),'not associated when do_FAR: ssi')
         _ASSERT(all(shape(ssi) == [ngptsw,ncol]),'mal-dimensioned: ssi')
         _ASSERT(associated(taucld_age),'not associated when do_FAR: taucld_age')
         _ASSERT(all(shape(taucld_age) == [ncol]),'mal-dimensioned: taucld_age')
         _ASSERT(associated(Rcldycol),'not associated when do_FAR: Rcldycol')
         _ASSERT(all(shape(Rcldycol) == [ncol]),'mal-dimensioned: Rcldycol')
         _ASSERT(associated(Rcldymc),'not associated when do_FAR: Rcldymc')
         _ASSERT(all(shape(Rcldymc) == [nlay,ngptsw,ncol]),'mal-dimensioned: Rcldymc')
         _ASSERT(associated(taucmc),'not associated when do_FAR: taucmc')
         _ASSERT(all(shape(taucmc) == [nlay,ngptsw,ncol]),'mal-dimensioned: taucmc')
         _ASSERT(associated(ssacmc),'not associated when do_FAR: ssacmc')
         _ASSERT(all(shape(ssacmc) == [nlay,ngptsw,ncol]),'mal-dimensioned: ssacmc')
         _ASSERT(associated(asmcmc),'not associated when do_FAR: asmcmc')
         _ASSERT(all(shape(asmcmc) == [nlay,ngptsw,ncol]),'mal-dimensioned: asmcmc')
         _ASSERT(associated(taormc),'not associated when do_FAR: taormc')
         _ASSERT(all(shape(taormc) == [nlay,ngptsw,ncol]),'mal-dimensioned: taormc')
      end if

      ! set column partition size pncol
      if (rpart > 0) then
         pncol = rpart
      else
         pncol = 2
      end if
      
      ! do partitions
      call rrtmg_sw_sub (MAPL, &
         pncol, ncol, nlay, &
         scon, adjes, coszen, isolvar, &
         play, plev, tlay, &
         h2ovmr, o3vmr, co2vmr, ch4vmr, o2vmr, &
         iceflgsw, liqflgsw, &
         cld, ciwp, clwp, rei, rel, &
         dyofyr, zm, alat, &
         iaer, tauaer, ssaaer, asmaer, &
         asdir, asdif, aldir, aldif, &
         cloudLM, cloudMH, normFlx, &
         clearCounts, swuflx, swdflx, swuflxc, swdflxc, &
         nirr, nirf, parr, parf, uvrr, uvrf, &
         tautp, tauhp, taump, taulp, &
         do_FAR, taumol_age, taumol_age_limit, &
         taur, taug, sflxzen, ssi, &
         taucld_age, taucld_age_limit, &
         Rcldycol, Rcldymc, taucmc, ssacmc, asmcmc, taormc, &
         bndscl, indsolvar, solcycfrac, &  ! optional inputs
         __RC__)
                                                      
      _RETURN(_SUCCESS)
   end subroutine rrtmg_sw                                                     


   subroutine rrtmg_sw_sub (MAPL, &
      pncol, gncol, nlay, &
      scon, adjes, gcoszen, isolvar, &
      gplay, gplev, gtlay, &
      gh2ovmr, go3vmr, gco2vmr, gch4vmr, go2vmr, &
      iceflgsw, liqflgsw, &
      gcld, gciwp, gclwp, grei, grel, &
      dyofyr, gzm, galat, &
      iaer, gtauaer, gssaaer, gasmaer, &
      gasdir, gasdif, galdir, galdif, &
      cloudLM, cloudMH, normFlx, &
      gclearCounts, gswuflx, gswdflx, gswuflxc, gswdflxc, &
      gnirr, gnirf, gparr, gparf, guvrr, guvrf, &
      gtautp, gtauhp, gtaump, gtaulp, &
      do_FAR, gtaumol_age, taumol_age_limit, &
      gtaur, gtaug, gsflxzen, gssi, &
      gtaucld_age, taucld_age_limit, &
      Rgcldycol, Rgcldymc, gtaucmc, gssacmc, gasmcmc, gtaormc, &
      bndscl, indsolvar, solcycfrac, &  ! optional inputs
      RC)


     ! ----- Modules -----
      use parrrsw, only : nbndsw, ngptsw, &
                          jpband, jpb1, jpb2, rrsw_scon
      use rrsw_con, only : grav, avogad
      use NRLSSI2, only : initialize_NRLSSI2, &
                          adjust_solcyc_amplitudes, &
                          interpolate_indices, &
                          Iint, Fint, Sint, &
                          Mg_avg, Mg_0, &
                          SB_avg, SB_0, &
                          isolvar_1_mean_svar_f, &
                          isolvar_1_mean_svar_s

      ! ----- Inputs -----
      ! (see rrtmg_sw() for more detailed comments)

      type(MAPL_MetaComp), intent(inout) :: MAPL

      ! dimensions
      integer, intent(in) :: pncol                     ! Nominal horiz cols in a partition
      integer, intent(in) :: gncol                     ! Global number of horizontal columns
      integer, intent(in) :: nlay                      ! Number of model layers

      ! orbit
      real, intent(in) :: scon                         ! Solar constant (W/m2)
      real, intent(in) :: adjes                        ! Flux adjustment for Earth/Sun distance
      real, intent(in) :: gcoszen (gncol)              ! Cosine of solar zenith angle

      ! solar variability
      integer, intent(in) :: isolvar
      real, intent(in), optional :: indsolvar (2)
      real, intent(in), optional :: bndscl (nbndsw)
      real, intent(in), optional :: solcycfrac

      ! profile
      real, intent(in) :: gplay   (gncol,nlay)         ! Layer pressures (hPa)
      real, intent(in) :: gplev   (gncol,nlay+1)       ! Interface pressures (hPa)
      real, intent(in) :: gtlay   (gncol,nlay)         ! Layer temperatures (K)

      ! gases
      real, intent(in) :: gh2ovmr (gncol,nlay)         ! H2O volume mixing ratio
      real, intent(in) :: go3vmr  (gncol,nlay)         ! O3 volume mixing ratio
      real, intent(in) :: gco2vmr (gncol,nlay)         ! CO2 volume mixing ratio
      real, intent(in) :: gch4vmr (gncol,nlay)         ! Methane volume mixing ratio
      real, intent(in) :: go2vmr  (gncol,nlay)         ! Oxygen volume mixing ratio

      ! cloud optics flags
      integer, intent(in) :: iceflgsw                  ! Flag for ice particle specifn
      integer, intent(in) :: liqflgsw                  ! Flag for liquid droplet specifn
      
      ! clouds
      real, intent(in) :: gcld    (gncol,nlay)         ! Cloud fraction
      real, intent(in) :: gciwp   (gncol,nlay)         ! In-cloud ice water path (g/m2)
      real, intent(in) :: gclwp   (gncol,nlay)         ! In-cloud liquid water path (g/m2)
      real, intent(in) :: grei    (gncol,nlay)         ! Cloud ice effective radius (um)
      real, intent(in) :: grel    (gncol,nlay)         ! Cloud drop effective radius (um)
                                                      
      ! cloud overlap
      integer, intent(in) :: dyofyr                    ! Day of the year
      real, intent(in) :: gzm     (gncol,nlay)         ! Heights of level midpoints
      real, intent(in) :: galat   (gncol)              ! Latitudes of columns [radians]
                                              
      ! aerosols (optical props, non-delta-scaled)
      integer, intent(in) :: iaer                      ! aerosol flag (0=off, 10=on)
      real, intent(in) :: gtauaer (gncol,nlay,nbndsw)  ! aer optical depth   (iaer=10 only)    
      real, intent(in) :: gssaaer (gncol,nlay,nbndsw)  ! aer single scat alb (iaer=10 only)    
      real, intent(in) :: gasmaer (gncol,nlay,nbndsw)  ! aer asymmetry param (iaer=10 only)    

      ! surface albedos
      real, intent(in) :: gasdir  (gncol)              ! UV/vis  surface albedo: direct rad
      real, intent(in) :: gasdif  (gncol)              ! UV/vis  surface albedo: diffuse rad
      real, intent(in) :: galdir  (gncol)              ! Near-IR surface albedo: direct rad
      real, intent(in) :: galdif  (gncol)              ! Near-IR surface albedo: diffuse rad

      ! super-layer cloud fraction boundaries 
      integer, intent(in) :: cloudLM                   ! Low-mid
      integer, intent(in) :: cloudMH                   ! Mid-high

      integer, intent(in) :: normFlx                   ! Normalize fluxes flag

      ! FAR controls
      logical, intent(in) :: do_FAR
      real, intent(in) :: taumol_age_limit
      real, intent(in) :: taucld_age_limit

      ! ----- Outputs -----

      ! subcolumn clear counts for Tot|High|Mid|Low super-layers
      integer, intent(out) :: gclearCounts(gncol,4)

      real, intent(out) :: gswuflx  (gncol,nlay+1)     !   All-sky SW up   flux (W/m2)
      real, intent(out) :: gswdflx  (gncol,nlay+1)     !   All-sky SW down flux (W/m2)
      real, intent(out) :: gswuflxc (gncol,nlay+1)     ! Clear-sky SW up   flux (W/m2)
      real, intent(out) :: gswdflxc (gncol,nlay+1)     ! Clear-sky SW down flux (W/m2)

      ! Output added for Land/Surface process (all-sky)
      real, intent(out) :: gnirr    (gncol)            ! Near-IR direct  down SW flux (w/m2)
      real, intent(out) :: gnirf    (gncol)            ! Near-IR diffuse down SW flux (w/m2)
      real, intent(out) :: gparr    (gncol)            ! Visible direct  down SW flux (w/m2)
      real, intent(out) :: gparf    (gncol)            ! Visible diffuse down SW flux (w/m2)
      real, intent(out) :: guvrr    (gncol)            ! UV      direct  down SW flux (w/m2)
      real, intent(out) :: guvrf    (gncol)            ! UV      diffuse down SW flux (w/m2)

      ! In-cloud PAR optical thickness for Tot|High|Mid|Low super-layers
      real, intent(out), dimension (gncol) :: gtautp, gtauhp, gtaump, gtaulp

      integer, intent(out), optional :: RC  ! return code

      ! ------- FAR InOuts -------
      ! if (.not.do_FAR) these can be unassociated pointers since not used

      real, intent(inout), dimension(:),     pointer :: gtaumol_age      !             (gncol) if (do_FAR)
      real, intent(inout), dimension(:,:,:), pointer :: gtaur, gtaug     ! (nlay,ngptsw,gncol) if (do_FAR)
      real, intent(inout), dimension(:,:),   pointer :: gsflxzen, gssi   ! (     ngptsw,gncol) if (do_FAR)
      real, intent(inout), dimension(:),     pointer :: gtaucld_age      !             (gncol) if (do_FAR)
      REAL, intent(inout), dimension(:),     pointer ::Rgcldycol         !             (gncol) if (do_FAR)
      REAL, intent(inout), dimension(:,:,:), pointer ::Rgcldymc          ! (nlay,ngptsw,gncol) if (do_FAR)
      real, intent(inout), dimension(:,:,:), pointer :: gtaucmc, gssacmc ! (nlay,ngptsw,gncol) if (do_FAR)
      real, intent(inout), dimension(:,:,:), pointer :: gasmcmc, gtaormc ! (nlay,ngptsw,gncol) if (do_FAR)

      ! ----- Locals -----

      ! Control
      real, parameter :: zepzen = 1.e-10  ! very small cossza
      integer :: ibnd, icol, ilay, ilev   ! various indices

      ! Atmosphere
      real :: coldry (nlay,pncol)     ! dry air column amount

      ! solar input
      real :: coszen (pncol)          ! Cosine of solar zenith angle
      real :: cossza (pncol)          ! Cosine of solar zenith angle
      real :: adjflux (jpband)        ! adjustment for curr Earth/Sun distance
      real :: gswdflx_at_top (gncol)  ! swdflx at TOA

      ! surface albedos
      real :: albdir (nbndsw,pncol)   ! surface albedo, direct
      real :: albdif (nbndsw,pncol)   ! surface albedo, diffuse
      
      ! Atmosphere/gases    
      ! ----------------

      ! general
      real :: play (nlay,  pncol)     ! Layer pressures (hPa)
      real :: plev (nlay+1,pncol)     ! Interface pressures (hPa)
      real :: tlay (nlay,  pncol)     ! Layer temperatures (K)

      ! gasesous absorbers
      real :: colh2o  (nlay,pncol)    ! column amount (h2o)
      real :: colco2  (nlay,pncol)    ! column amount (co2)
      real :: colo3   (nlay,pncol)    ! column amount (o3)
      real :: colch4  (nlay,pncol)    ! column amount (ch4)
      real :: colo2   (nlay,pncol)    ! column amount (o2)

      ! Atmosphere/clouds - cldprop
      ! ---------------------------

      integer :: ncbands                    ! num of cloud spectral bands

      real :: cld  (nlay,pncol)             ! Cloud fraction
      real :: ciwp (nlay,pncol)             ! In-cloud ice water path [g/m2]
      real :: clwp (nlay,pncol)             ! In-cloud liq water path [g/m2]
      real :: rei  (nlay,pncol)             ! Cloud ice effective radius [um]
      real :: rel  (nlay,pncol)             ! Cloud drop effective radius [um]
      
      real :: alat      (pncol)             ! latitude for cloud overlap
      real :: zm   (nlay,pncol)		    ! mid-layer hgt for cld overlap [m]
                                                      
      logical :: cldymc (nlay,ngptsw,pncol) ! cloud or not? [mcica]
      integer :: clearCounts (4,pncol)      ! for super-layer cld fractions

      real :: taucmc (nlay,ngptsw,pncol)    ! in-cloud optical depth [mcica]
      real :: taormc (nlay,ngptsw,pncol)    ! unscaled in-cloud optl depth [mcica]
      real :: ssacmc (nlay,ngptsw,pncol)    ! in-cloud single scat albedo [mcica]
      real :: asmcmc (nlay,ngptsw,pncol)    ! in-cloud asymmetry param [mcica]
      
      ! Atmosphere/clouds/aerosol - spcvrt,spcvmc
      ! -----------------------------------------

      real :: taua (nlay,nbndsw,pncol)
      real :: asya (nlay,nbndsw,pncol)
      real :: omga (nlay,nbndsw,pncol)

      ! SW flux temporaries [W/m2]
      real :: zbbfu    (nlay+1,pncol)  ! all-SW  up           all-sky
      real :: zbbfd    (nlay+1,pncol)  ! all-SW  down         all-sky
      real :: zbbcu    (nlay+1,pncol)  ! all-SW  up          clear-sky
      real :: zbbcd    (nlay+1,pncol)  ! all-SW  down        clear-sky
      real :: zbbfddir (nlay+1,pncol)  ! all-SW  down direct   all-sky
      real :: zbbcddir (nlay+1,pncol)  ! all-SW  down direct clear-sky

      real, dimension (pncol) :: &
         nirr, nirf, parr, parf, uvrr, uvrf

      ! in-cloud PAR optical thicknesses
      real, dimension (pncol) :: tautp, tauhp, taump, taulp
      
      ! FAR taumol partitioned fields
      real, dimension(pncol) :: tmage
      real, dimension(nlay,ngptsw,pncol) :: taur, taug
      real, dimension(ngptsw,pncol) :: sflxzen, ssi

      ! FAR taucld partitioned fields
      logical, dimension(pncol) :: tcrecalc

      ! Solar variability multipliers
      ! -----------------------------
      real :: svar_f               ! facular multiplier
      real :: svar_s               ! sunspot multiplier
      real :: svar_i               ! baseline irradiance multiplier
      real :: svar_f_bnd (jpband)  ! facular multiplier (by band)
      real :: svar_s_bnd (jpband)  ! sunspot multiplier (by band)
      real :: svar_i_bnd (jpband)  ! baseline irradiance multiplier (by band)

      integer :: ncol_clr, ncol_cld
      integer :: npart_clr, npart_cld, npart
      integer, dimension (gncol) :: gicol_clr, gicol_cld

      real, parameter :: amd = 28.9660     ! Effective molecular weight of dry air (g/mol)
      real, parameter :: amw = 18.0160     ! Molecular weight of water vapor (g/mol)

! Set molecular weight ratios (for converting mmr to vmr), e.g. h2ovmr = h2ommr * amdw

      real, parameter :: amdw  = 1.607793  ! Molecular weight of dry air / water vapor
      real, parameter :: amdc  = 0.658114  ! Molecular weight of dry air / carbon dioxide
      real, parameter :: amdo  = 0.603428  ! Molecular weight of dry air / ozone
      real, parameter :: amdm  = 1.805423  ! Molecular weight of dry air / methane
      real, parameter :: amdn  = 0.658090  ! Molecular weight of dry air / nitrous oxide
      real, parameter :: amdo2 = 0.905140  ! Molecular weight of dry air / oxygen

      integer :: n, imol, gicol            ! Loop indices
      real :: adjflx                       ! flux adjustment for Earth/Sun distance
      
      integer :: ipart, col_last, cols, cole, cc

      ! ncol is the actual number of gridcols in a partition, cf. pncol,
      ! the maximum number. May have ncol < pncol on final partition.
      integer :: ncol
      integer, allocatable :: idx(:)

      ! other solar variability locals
      ! ------------------------------
      real :: solvar (jpband)              ! solar constant scaling factor by band
      real :: indsolvar_scl (2)            ! Adjusted facular and sunspot amplitude 
                                           !   scale factors (isolvar=1)
      real :: indsolvar_ndx (2)            ! Facular and sunspot indices (isolvar=2)

      real :: solcycfr, Mg_now, SB_now
      real :: scon_int, svar_r
 
      ! FAR locals
      logical :: gtaucld_recalc (gncol)
      integer :: nrc
      integer :: irc (pncol)
      real :: alat_rc (pncol)
      integer :: clearCounts_rc(4,pncol)
      logical :: cldymc_rc (nlay,ngptsw,pncol)
      real    :: ciwpmc_rc (nlay,ngptsw,pncol)  ! in-cloud ice water path [mcica] [g/m2]
      real    :: clwpmc_rc (nlay,ngptsw,pncol)  ! in-cloud liq water path [mcica] [g/m2]
      real, dimension(nlay,pncol) :: &
         zm_rc, play_rc, cld_rc, ciwp_rc, clwp_rc, rei_rc, rel_rc
      real, dimension(nlay,ngptsw,pncol) :: &
         taucmc_rc, ssacmc_rc, asmcmc_rc, taormc_rc

      integer :: STATUS  ! for MAPL error reporting

      ! Having no work to do at all would mess up timings
      ! This passes because load balancer is multi-pass.
      _ASSERT(gncol > 0, 'no columns on processor!')

      ! Initializations
      ! ---------------

      ! solar variability: default values
      solvar(:) = 1.
      adjflux(:) = 1.
      svar_f = 1.
      svar_s = 1. 
      svar_i = 1. 
      svar_f_bnd(:) = 1. 
      svar_s_bnd(:) = 1. 
      svar_i_bnd(:) = 1. 

      ! isolvar == 1 specifies the position in AvgCyc11 through solcycfrac
      ! and allows scaling of solar cycle amplitudes as described in notes.
      ! ------------------------------------------------------------------

      if (isolvar .eq. 1) then 

         ! require solcycfrac present, else what's the point of using isolvar=1 ?
         if (.not.present(solcycfrac)) then
            _FAIL('isolvar == 1 requires solcycfrac present!')
         end if
         solcycfr = solcycfrac

         ! No amplitude scaling unless indsolvar is present. 
         indsolvar_scl(1:2) = 1.

         if (present(indsolvar)) then 

            ! Adjust amplitude scaling of mean solar cycle to be unity at
            ! solar minimum (solcycfrac_min), to be the requested indsolvar
            ! at solar maximum (solcycfrac_max), and to vary linearly with
            ! solcycfr between those values.

            if (indsolvar(1).ne.1. .or. indsolvar(2).ne.1.) &
               call adjust_solcyc_amplitudes(solcycfr, indsolvar, indsolvar_scl)

         endif

      endif

      ! isolvar == 2 allows direct specification of Mg and SB via indsolvar
      ! -------------------------------------------------------------------
      
      if (isolvar .eq. 2) then 

         ! default to mean indices
         indsolvar_ndx(1) = Mg_avg
         indsolvar_ndx(2) = SB_avg

         ! update to specified indices if provided
         if (present(indsolvar)) then 
            indsolvar_ndx(1) = indsolvar(1)
            indsolvar_ndx(2) = indsolvar(2)
         endif

      endif

      ! pre-calculated constants (will only do calcs once internally)
      ! -------------------------------------------------------------
      call initialize_NRLSSI2 (isolvar, indsolvar)

      ! Set flux adjustment for current Earth/Sun distance (two options)
      ! ----------------------------------------------------------------
      ! (Set adjflx to 1. to use constant Earth/Sun distance of 1 AU). 

      ! 1) Provided by GCM via ADJES (from MAPL sun factor DIST ~ 1/r^2)
      adjflx = adjes

      ! 2) Calc Earth/Sun dist adj from DYOFYR, the cumulative day of year
      ! (Turned off but DYOFYR used by MCICA exponential cloud overlap).

      ! if (dyofyr .gt. 0) then
      !    adjflx = earth_sun(dyofyr)
      ! endif

      ! --------------------------------------------------------
      ! Apply selected solar variability option based on ISOLVAR
      ! and input solar constant SCON.
      ! --------------------------------------------------------

      if (scon == 0.) then 

         ! For scon = 0, use internally defined solar constant, which is
         ! 1368.22 Wm-2 (for ISOLVAR=-1) and 1360.85 Wm-2 (For ISOLVAR=0,3;
         ! Options ISOLVAR=1,2 model sol cyc varations from 1360.85 Wm-2).

         if (isolvar .eq. -1) then

            ! Constant sun (Kurucz)
            ! Apply optional scaling by band if bndscl present.

            solvar(jpb1:jpb2) = 1.
            if (present(bndscl)) solvar(jpb1:jpb2) = bndscl(:)

         elseif (isolvar .eq. 0) then

            ! Constant sun (NRLSSI2 model)
            ! Quiet sun, facular, and sunspot terms averaged over AvgCyc11.

            svar_f = 1.
            svar_s = 1.
            svar_i = 1.

         elseif (isolvar .eq. 1) then

            ! Apply NRLSSI2 solar irradiance model at a specified solcycfr
            ! within AvgCyc11, with the additional amplitude scalings in 
            ! indsolvar_scl.

            ! interpolate mean solar cycle to solcycfr
            call interpolate_indices (solcycfr, Mg_now, SB_now)

            ! Apply linear index-to-flux-multiplier-svar relationship
            ! with the additional indsolvar_scl scaling.
            svar_f = indsolvar_scl(1) * (Mg_now - Mg_0) / (Mg_avg - Mg_0)
            svar_s = indsolvar_scl(2) * (SB_now - SB_0) / (SB_avg - SB_0)
            svar_i = 1.

         elseif (isolvar .eq. 2) then

            ! Specified solar cycle with solar variability based on NRLSSI2 model.
            ! Facular and sunspot index terms input directly.

            svar_f = (indsolvar_ndx(1) - Mg_0) / (Mg_avg - Mg_0)
            svar_s = (indsolvar_ndx(2) - SB_0) / (SB_avg - SB_0)
            svar_i = 1.

         elseif (isolvar .eq. 3) then

            ! Constant sun (NRLSSI2 model).
            ! Averaged facular, sunspot and quiet sun terms from AvgCyc11.
            ! Apply optional scaling by band if bndscl present.

            solvar(jpb1:jpb2) = 1.
            if (present(bndscl)) solvar(jpb1:jpb2) = bndscl(:)
            do ibnd = jpb1,jpb2
               svar_f_bnd(ibnd) = solvar(ibnd)
               svar_s_bnd(ibnd) = solvar(ibnd)
               svar_i_bnd(ibnd) = solvar(ibnd)
            enddo

         else
            _FAIL('invalid isolvar')
         endif 

      elseif (scon > 0.) then 

         ! Scale from internal to externally specified SCON.

         if (isolvar .eq. -1) then

            ! Constant sun (Kurucz)
            ! Scale from internal to requested solar constant.
            ! Apply optional scaling by band if bndscl present.

            solvar(jpb1:jpb2) = scon / rrsw_scon 
            if (present(bndscl)) &
               solvar(jpb1:jpb2) = solvar(jpb1:jpb2) * bndscl(:)

         elseif (isolvar .eq. 0) then

            ! Constant sun (NRLSSI2 model)
            ! Quiet sun, facular, and sunspot terms averaged over AvgCyc11.
            ! Scale from internal to requested solar constant. 

            scon_int = Fint + Sint + Iint
            svar_r = scon / scon_int
            svar_f = svar_r
            svar_s = svar_r
            svar_i = svar_r

         elseif (isolvar .eq. 1) then

            ! Apply NRLSSI2 solar irradiance model at a specified solcycfr
            ! within AvgCyc11, with the additional amplitude scalings in
            ! indsolvar_scl. Scale from the internal to the requested solar
            ! constant, which is treated as a required *cycle average*.

            ! interpolate mean solar cycle to solcycfr
            call interpolate_indices (solcycfr, Mg_now, SB_now)

            ! Apply linear index-to-flux-multiplier-svar relationship
            ! with the additional indsolvar_scl scaling. Select a constant
            ! svar_i such that chosen scon is the <cycle average>.
            ! scon = svar_i * Iint + <svar_f> * Fint + <svar_s> * Sint >
            ! => svar_i = [scon - (<svar_f> * Fint + <svar_s> * Sint)] / Iint

            svar_f = indsolvar_scl(1) * (Mg_now - Mg_0) / (Mg_avg - Mg_0)
            svar_s = indsolvar_scl(2) * (SB_now - SB_0) / (SB_avg - SB_0)
            svar_i = (scon - (isolvar_1_mean_svar_f * Fint + &
                              isolvar_1_mean_svar_s * Sint)) / Iint

         elseif (isolvar .eq. 2) then

            ! Specified solar cycle with solar variability based on NRLSSI2 model.
            ! Facular and sunspot index terms input directly. Scale from internal
            ! to requested solar constant by setting svar_i so that
            !   svar_i * Iint + svar_f * Fint + svar_s * Sint = scon.
            ! So, scon is honored at EACH time, because it too is assumed to be
            ! specified from time-varying data.

            svar_f = (indsolvar_ndx(1) - Mg_0) / (Mg_avg - Mg_0)
            svar_s = (indsolvar_ndx(2) - SB_0) / (SB_avg - SB_0)
            svar_i = (scon - (svar_f * Fint + svar_s * Sint)) / Iint 

         elseif (isolvar .eq. 3) then

            ! Constant sun (NRLSSI2 model).
            ! Averaged facular, sunspot and quiet sun terms from AvgCyc11.
            ! Scale from internal to requested solar constant.
            ! Apply optional scaling by band if bndscl present.

            scon_int = Fint + Sint + Iint
            solvar(jpb1:jpb2) = scon / scon_int
            if (present(bndscl)) solvar(jpb1:jpb2) = solvar(jpb1:jpb2) * bndscl(:)
            do ibnd = jpb1,jpb2
               svar_f_bnd(ibnd) = solvar(ibnd)
               svar_s_bnd(ibnd) = solvar(ibnd)
               svar_i_bnd(ibnd) = solvar(ibnd)
            enddo

         else
            _FAIL('invalid isolvar')
         endif 

      else
         _FAIL('scon cannot be negative!')
      endif

      ! Earth-Sun distance adjustment
      adjflux(jpb1:jpb2) = adjflx

      ! Combine with solar constant scaling for Kurucz
      ! (done separately via svar_ for NRLSSI2)
      if (isolvar < 0) then
         adjflux(jpb1:jpb2) = adjflux(jpb1:jpb2) * solvar(jpb1:jpb2)
      endif

      ! Which columns to recalculate cloud optical properties?
      ! Needed, if FAR, before clear/cloudy separation below.
      if (.not.do_FAR) then
         ! all of them
         gtaucld_recalc = .true.
      else
         ! FAR: asynchronous recalculation of uninitialized or old values ...
         gtaucld_recalc = (gtaucld_age < 0. .or. gtaucld_age > taucld_age_limit)
         ! Set soon-to-be recalculated values to brand new.
         where (gtaucld_recalc) gtaucld_age = 0.
      endif

      ! Build profile separation based on cloudiness, i.e., count and index
      ! clear/cloudy gridcolumns. The separation is based on whether the grid-
      ! column has cloud fraction in any layer (or not). This is based on the
      ! gridcolumn state, not on the later generated McICA subcolumn ensemble.
      ! If all layers have zero cld fraction then so will all McICA subcolumns,
      ! but the converse in not true ... can easily get a clear subcolumn for a
      ! cloudy gridcolumn. So, the gicol_clr can be assumed to yield all clear
      ! subcolumns, while the gicol_cld will yield both clear and cloudy sub-
      ! columns. 
      ! Also, if FAR, recalculate column cloudiness as necessary.
      ncol_clr = 0
      ncol_cld = 0
      if (.not.do_FAR) then
         do gicol = 1,gncol
            if (any(gcld(gicol,:) > 0)) then
               ncol_cld = ncol_cld + 1
               gicol_cld(ncol_cld) = gicol
            else
               ncol_clr = ncol_clr + 1
               gicol_clr(ncol_clr) = gicol
            end if
         end do
      else  ! FAR
         do gicol = 1,gncol
            if (gtaucld_recalc(gicol)) Rgcldycol(gicol) = merge(1., 0., any(gcld(gicol,:) > 0))
            if (Rgcldycol(gicol).ne.0.) then
               ncol_cld = ncol_cld + 1
               gicol_cld(ncol_cld) = gicol
            else
               ncol_clr = ncol_clr + 1
               gicol_clr(ncol_clr) = gicol
            end if
         end do
      end if

      ! num of length pncol partitions needed for clear/cloudy profiles
      npart_clr = ceiling( real(ncol_clr) / real(pncol) )
      npart_cld = ceiling( real(ncol_cld) / real(pncol) )

      ! "zero" aerosols if dont want them included
      if (iaer /= 10) then
         taua = 0.
         asya = 0.
         omga = 1.
      end if

      ! partitioning over clear (cc=1) and cloudy (cc=2) gridcolumns
      ! ------------------------------------------------------------

      do cc = 1,2  ! outer loop over clear then cloudy gridcolumns

         if (cc == 1) then 
            ! clear
            npart = npart_clr
            col_last = ncol_clr
         else
            ! cloudy
            npart = npart_cld
            col_last = ncol_cld
         end if

         ! loop over partitions
         do ipart = 0,npart-1

            call MAPL_TimerOn (MAPL,"RRTMG_PART",__RC__)

            ! partition dimensions
            cols = ipart * pncol + 1
            cole = (ipart + 1) * pncol
            if (cole > col_last) cole = col_last
            ncol = cole - cols + 1
            allocate(idx(ncol),__STAT__)

            ! copy inputs into partition
            ! --------------------------

            if (cc == 1) then    

               ! -----------------
               ! Clear gridcolumns
               ! -----------------

               do icol = 1,ncol
                  gicol = gicol_clr(icol + cols - 1)

                  ! assign surface albedos to bands

                  ! near IR bands 14=nbndsw and 1-8
                  ! 820-12850 cm-1, 0.778-12.2 um
                  do ibnd=1,8
                     albdir(ibnd,icol) = galdir(gicol)
                     albdif(ibnd,icol) = galdif(gicol)
                  enddo
                  albdir(nbndsw,icol) = galdir(gicol)
                  albdif(nbndsw,icol) = galdif(gicol)

                  ! UV/Vis bands 10-13
                  ! 16000-50000 cm-1, 0.200-0.625 um
                  do ibnd=10,13
                     albdir(ibnd,icol) = gasdir(gicol)
                     albdif(ibnd,icol) = gasdif(gicol)
                  enddo

                  ! Transition band 9
                  ! 12850-16000 cm-1, 0.625-0.778 um
                  ! Take average, dmlee
                  albdir(9,icol) = (gasdir(gicol)+galdir(gicol))/2.
                  albdif(9,icol) = (gasdif(gicol)+galdif(gicol))/2.

               enddo

               ! copy in partition (general)
               do icol = 1,ncol
                  gicol = gicol_clr(icol + cols - 1)
    
                  play(:,icol) = gplay(gicol,1:nlay)
                  plev(:,icol) = gplev(gicol,1:nlay+1)
                  tlay(:,icol) = gtlay(gicol,1:nlay)
                  coszen(icol) = gcoszen(gicol)

               enddo

               ! copy in partition (aerosols)
               if (iaer == 10) then
                  do icol = 1,ncol
                     gicol = gicol_clr(icol + cols - 1)
                     do ibnd = 1,nbndsw
                        taua(1:nlay,ibnd,icol) = gtauaer(gicol,1:nlay,ibnd)
                        asya(1:nlay,ibnd,icol) = gasmaer(gicol,1:nlay,ibnd)
                        omga(1:nlay,ibnd,icol) = gssaaer(gicol,1:nlay,ibnd)
                     enddo
                  enddo
               endif   

               ! copy in partition (gases)
               do icol = 1,ncol
                  gicol = gicol_clr(icol + cols - 1)
                  colh2o(:,icol) = gh2ovmr(gicol,1:nlay)
                  colco2(:,icol) = gco2vmr(gicol,1:nlay)
                  colo3 (:,icol) = go3vmr (gicol,1:nlay)
                  colch4(:,icol) = gch4vmr(gicol,1:nlay)
                  colo2 (:,icol) = go2vmr (gicol,1:nlay)   
                end do

!pmn:??speed investigate efficiency in other copy-in/out sections

               ! copy in FAR taumol InOuts:
               if (do_FAR) then
                  idx = gicol_clr(cols:cole)
                  tmage     (1:ncol) =  gtaumol_age(idx)
                  sflxzen (:,1:ncol) =  gsflxzen (:,idx)
                  ssi     (:,1:ncol) =  gssi     (:,idx)
                  taur  (:,:,1:ncol) =  gtaur  (:,:,idx)
                  taug  (:,:,1:ncol) =  gtaug  (:,:,idx)
                  tcrecalc  (1:ncol) =  gtaucld_recalc(idx)
                  cldymc(:,:,1:ncol) =(Rgcldymc(:,:,idx).ne.0.)
                  taucmc(:,:,1:ncol) =  gtaucmc(:,:,idx)
                  ssacmc(:,:,1:ncol) =  gssacmc(:,:,idx)
                  asmcmc(:,:,1:ncol) =  gasmcmc(:,:,idx)
                  taormc(:,:,1:ncol) =  gtaormc(:,:,idx)
               end if

            else

               ! ------------------
               ! Cloudy gridcolumns
               ! ------------------
          
               do icol = 1,ncol
                  gicol = gicol_cld(icol + cols - 1)
     
                  ! assign surface albedos to bands

                  ! near IR bands 14=nbndsw and 1-8
                  ! 820-12850 cm-1, 0.778-12.2 um
                  do ibnd=1,8
                     albdir(ibnd,icol) = galdir(gicol)
                     albdif(ibnd,icol) = galdif(gicol)
                  enddo
                  albdir(nbndsw,icol) = galdir(gicol)
                  albdif(nbndsw,icol) = galdif(gicol)

                  ! UV/Vis bands 10-13
                  ! 16000-50000 cm-1, 0.200-0.625 um
                  do ibnd=10,13
                     albdir(ibnd,icol) = gasdir(gicol)
                     albdif(ibnd,icol) = gasdif(gicol)
                  enddo

                  ! Transition band 9
                  ! 12850-16000 cm-1, 0.625-0.778 um
                  ! Take average, dmlee
                  albdir(9,icol) = (gasdir(gicol)+galdir(gicol))/2.
                  albdif(9,icol) = (gasdif(gicol)+galdif(gicol))/2.

               enddo
          
               ! copy in partition (general and cloud physical props)
               do icol = 1,ncol
                  gicol = gicol_cld(icol + cols - 1)
     
                  play(:,icol) = gplay(gicol,1:nlay)
                  plev(:,icol) = gplev(gicol,1:nlay+1)
                  tlay(:,icol) = gtlay(gicol,1:nlay)
                  cld (:,icol) = gcld (gicol,1:nlay)
                  ciwp(:,icol) = gciwp(gicol,1:nlay)
                  clwp(:,icol) = gclwp(gicol,1:nlay)
                  rei (:,icol) = grei (gicol,1:nlay) 
                  rel (:,icol) = grel (gicol,1:nlay)
                  zm  (:,icol) = gzm  (gicol,1:nlay)
                  alat  (icol) = galat  (gicol)
                  coszen(icol) = gcoszen(gicol)
               enddo

               ! copy in partition (aerosols)
               if (iaer == 10) then
                  do icol = 1,ncol
                     gicol = gicol_cld(icol + cols - 1)
                     do ibnd = 1,nbndsw
                        taua(1:nlay,ibnd,icol) = gtauaer(gicol,1:nlay,ibnd)
                        asya(1:nlay,ibnd,icol) = gasmaer(gicol,1:nlay,ibnd)
                        omga(1:nlay,ibnd,icol) = gssaaer(gicol,1:nlay,ibnd)
                     end do
                  end do
               endif

               ! copy in partition (gases)
               do icol = 1,ncol
                  gicol = gicol_cld(icol + cols - 1)
                  colh2o(:,icol) = gh2ovmr(gicol,1:nlay)
                  colco2(:,icol) = gco2vmr(gicol,1:nlay)
                  colo3 (:,icol) = go3vmr (gicol,1:nlay)
                  colch4(:,icol) = gch4vmr(gicol,1:nlay)
                  colo2 (:,icol) = go2vmr (gicol,1:nlay)  
               enddo

               ! copy in FAR taumol InOuts
               if (do_FAR) then
                  idx = gicol_cld(cols:cole)
                  tmage     (1:ncol) =  gtaumol_age(idx)
                  sflxzen (:,1:ncol) =  gsflxzen (:,idx)
                  ssi     (:,1:ncol) =  gssi     (:,idx)
                  taur  (:,:,1:ncol) =  gtaur  (:,:,idx)
                  taug  (:,:,1:ncol) =  gtaug  (:,:,idx)
                  tcrecalc  (1:ncol) =  gtaucld_recalc(idx)
                  cldymc(:,:,1:ncol) =(Rgcldymc(:,:,idx).ne.0.)
                  taucmc(:,:,1:ncol) =  gtaucmc(:,:,idx)
                  ssacmc(:,:,1:ncol) =  gssacmc(:,:,idx)
                  asmcmc(:,:,1:ncol) =  gasmcmc(:,:,idx)
                  taormc(:,:,1:ncol) =  gtaormc(:,:,idx)
               end if

            end if  ! clear or cloudy gridcolumns

            call MAPL_TimerOff(MAPL,"RRTMG_PART",__RC__)

            ! limit tiny cosine zenith angles
            do icol = 1,ncol
               cossza(icol) = max(zepzen,coszen(icol))
            enddo

            ! evaluate dry air molecules/cm^2
            ! (see details in rrtmg_lw_rad())
            do icol = 1,ncol
               do ilay = 1,nlay
                  coldry(ilay,icol) = (plev(ilay,icol)-plev(ilay+1,icol)) * 1.e3 * avogad / &
                     (1.e2 * grav * ((1.-colh2o(ilay,icol)) * amd + colh2o(ilay,icol) * amw) * &
                     (1. + colh2o(ilay,icol)))
               enddo
            enddo

            ! gases also to molecules/cm^2
            do icol = 1,ncol
               do ilay = 1,nlay
                  colh2o(ilay,icol) = coldry(ilay,icol) * colh2o(ilay,icol)
                  colco2(ilay,icol) = coldry(ilay,icol) * colco2(ilay,icol)
                  colo3 (ilay,icol) = coldry(ilay,icol) * colo3 (ilay,icol)
                  colch4(ilay,icol) = coldry(ilay,icol) * colch4(ilay,icol)
                  colo2 (ilay,icol) = coldry(ilay,icol) * colo2 (ilay,icol)
               end do
            end do

            ! cloudy gridcolumns
            if (cc == 2) then

               ! which columns to recalculate?
               if (.not.do_FAR) then

                  ! all of them
                  nrc = ncol
                  irc(1:ncol) = [1:ncol]

               else  ! FAR: asynchronous recalculation of uninitialized or old values ...

                  ! Get number of recalculated columns and their indicies irc.
                  nrc = 0
                  do icol = 1,ncol
                     if (tcrecalc(icol)) then
                        nrc = nrc + 1
                        irc(nrc) = icol
                     end if
                  end do

               endif

               ! recalculate as needed
               if (nrc > 0) then

                  ! McICA subcolumn generation
                  call MAPL_TimerOn (MAPL,"RRTMG_CLDSGEN",__RC__)
                  ! copy-in inputs for recalculated columns
                  do n = 1,nrc
                     icol = irc(n)
                     alat_rc(  n) = alat(  icol)
                       zm_rc(:,n) =   zm(:,icol)
                     play_rc(:,n) = play(:,icol)
                      cld_rc(:,n) =  cld(:,icol)
                     ciwp_rc(:,n) = ciwp(:,icol)
                     clwp_rc(:,n) = clwp(:,icol)
                  end do
                  call generate_stochastic_clouds( &
                     pncol, nrc, ngptsw, nlay, &
                     zm_rc, alat_rc, dyofyr, &
                     play_rc, cld_rc, ciwp_rc, clwp_rc, 1.e-20, &
                     cldymc_rc, ciwpmc_rc, clwpmc_rc, &
                     seed_order=[4,3,2,1]) 

! pmn: idea ... save space later by cldymc = ciwp > 0. .or. clwp > 0
! and then transfer to taucmc > 0.

                  ! for super-layer cloud fractions
                  call clearCounts_threeBand( &
                     pncol, nrc, ngptsw, nlay, cloudLM, cloudMH, cldymc_rc, &
                     clearCounts_rc)
                  ! copy-out recalculated values
                  do n = 1,nrc
                     icol = irc(n)
                     cldymc   (:,:,icol) =      cldymc_rc(:,:,n)
                     clearCounts(:,icol) = clearCounts_rc(  :,n)
                  end do
                  call MAPL_TimerOff(MAPL,"RRTMG_CLDSGEN",__RC__)

                  ! cloud optical property generation
                  call MAPL_TimerOn (MAPL,"RRTMG_CLDPRMC",__RC__)
                  ! copy-in inputs for recalculated columns
                  do n = 1,nrc
                     icol = irc(n)
                     rei_rc(:,n) = rei(:,icol)
                     rel_rc(:,n) = rel(:,icol)
                  end do
                  call cldprmc_sw( &
                     pncol, nrc, nlay, iceflgsw, liqflgsw,  &
                     cldymc_rc, ciwpmc_rc, clwpmc_rc, rei_rc, rel_rc, &
                     taormc_rc, taucmc_rc, ssacmc_rc, asmcmc_rc)
                  ! copy-out recalculated values
                  do n = 1,nrc
                     icol = irc(n)
                     taucmc(:,:,icol) = taucmc_rc(:,:,n)
                     ssacmc(:,:,icol) = ssacmc_rc(:,:,n)
                     asmcmc(:,:,icol) = asmcmc_rc(:,:,n)
                     taormc(:,:,icol) = taormc_rc(:,:,n)
                  end do
                  call MAPL_TimerOff(MAPL,"RRTMG_CLDPRMC",__RC__)

               else  ! nrc == 0
                  call MAPL_TimerOn (MAPL,"RRTMG_CLDSGEN",__RC__)
                  call MAPL_TimerOff(MAPL,"RRTMG_CLDSGEN",__RC__)
                  call MAPL_TimerOn (MAPL,"RRTMG_CLDPRMC",__RC__)
                  call MAPL_TimerOff(MAPL,"RRTMG_CLDPRMC",__RC__)
               endif

            else  ! cc == 1 (clear gridcolumns)

               call MAPL_TimerOn (MAPL,"RRTMG_CLDSGEN",__RC__)
               do icol = 1,ncol
                  cldymc(:,:,icol) = .false.
               end do
               call MAPL_TimerOff(MAPL,"RRTMG_CLDSGEN",__RC__)
               call MAPL_TimerOn (MAPL,"RRTMG_CLDPRMC",__RC__)
               do icol = 1,ncol
                  taucmc(:,:,icol) = 0.
                  ssacmc(:,:,icol) = 1. 
                  asmcmc(:,:,icol) = 0.
                  taormc(:,:,icol) = 0.
               end do
               call MAPL_TimerOff(MAPL,"RRTMG_CLDPRMC",__RC__)
            end if

            ! compute sw radiative fluxes
            call spcvmc_sw (MAPL, &
               cc, pncol, ncol, nlay, &
               play, tlay, coldry, &
               albdif, albdir, &
               cldymc, taucmc, asmcmc, ssacmc, taormc, &
               taua, asya, omga, cossza, adjflux, &
               isolvar, svar_f, svar_s, svar_i, &
               svar_f_bnd, svar_s_bnd, svar_i_bnd, &
               colch4, colco2, colh2o, colo2, colo3, &
               cloudLM, cloudMH, & 
               zbbfd, zbbfu, zbbcd, zbbcu, zbbfddir, zbbcddir, &
               nirr, nirf, parr, parf, uvrr, uvrf, &
               tautp, tauhp, taump, taulp, &
               do_FAR, tmage, taumol_age_limit, &
               taur, taug, sflxzen, ssi, &
               __RC__)

            ! Copy out up and down, clear- and all-sky fluxes to output arrays.
            ! Vertical indexing goes from bottom to top; reverse here for GCM if necessary.

            call MAPL_TimerOn (MAPL,"RRTMG_PART",__RC__)

            if (cc == 1) then  ! clear gridcolumns

               do icol = 1,ncol
                  gicol = gicol_clr(icol + cols - 1)
        
                  ! super-layer clear counts
                  do n = 1,4
                     gclearCounts(gicol,n) = ngptsw
                  end do

                  ! up and down fluxes
                  do ilev = 1,nlay+1
                     gswuflxc(gicol,ilev) = zbbcu(ilev,icol) 
                     gswdflxc(gicol,ilev) = zbbcd(ilev,icol) 
                     gswuflx (gicol,ilev) = zbbfu(ilev,icol) 
                     gswdflx (gicol,ilev) = zbbfd(ilev,icol) 
                  enddo

                  ! super-layer optical thicknesses
                  gtautp(gicol) = 0.
                  gtauhp(gicol) = 0.
                  gtaump(gicol) = 0.
                  gtaulp(gicol) = 0.

               enddo

               ! surface broadband fluxes
               do icol = 1,ncol
                  gicol = gicol_clr(icol + cols - 1)
                  gnirr(gicol) = nirr(icol)
                  gnirf(gicol) = nirf(icol) - nirr(icol)
                  gparr(gicol) = parr(icol)
                  gparf(gicol) = parf(icol) - parr(icol)
                  guvrr(gicol) = uvrr(icol)
                  guvrf(gicol) = uvrf(icol) - uvrr(icol)
               end do

               ! copy out FAR taumol InOuts
               if (do_FAR) then
                  idx = gicol_clr(cols:cole)
                  gtaumol_age(idx) = tmage     (1:ncol)
                  gsflxzen (:,idx) = sflxzen (:,1:ncol)
                  gssi     (:,idx) = ssi     (:,1:ncol)
                  gtaur  (:,:,idx) = taur  (:,:,1:ncol)
                  gtaug  (:,:,idx) = taug  (:,:,1:ncol)
                 Rgcldymc(:,:,idx) = merge(1.,0.,cldymc(:,:,1:ncol))
                  gtaucmc(:,:,idx) = taucmc(:,:,1:ncol)
                  gssacmc(:,:,idx) = ssacmc(:,:,1:ncol)
                  gasmcmc(:,:,idx) = asmcmc(:,:,1:ncol)
                  gtaormc(:,:,idx) = taormc(:,:,1:ncol)
               end if

            else ! cloudy columns

               do icol = 1,ncol
                  gicol = gicol_cld(icol + cols - 1)
                  do n = 1,4
                     gclearCounts(gicol,n) = clearCounts(n,icol)
                  end do
                  do ilev = 1,nlay+1
                     gswuflxc(gicol,ilev) = zbbcu(ilev,icol) 
                     gswdflxc(gicol,ilev) = zbbcd(ilev,icol) 
                     gswuflx (gicol,ilev) = zbbfu(ilev,icol) 
                     gswdflx (gicol,ilev) = zbbfd(ilev,icol) 
                  enddo
                  gtautp(gicol) = tautp(icol)
                  gtauhp(gicol) = tauhp(icol)
                  gtaump(gicol) = taump(icol)
                  gtaulp(gicol) = taulp(icol)
               enddo

               do icol = 1,ncol
                  gicol = gicol_cld(icol + cols - 1)
                  gnirr(gicol) = nirr(icol)
                  gnirf(gicol) = nirf(icol) - nirr(icol)
                  gparr(gicol) = parr(icol)
                  gparf(gicol) = parf(icol) - parr(icol)
                  guvrr(gicol) = uvrr(icol)
                  guvrf(gicol) = uvrf(icol) - uvrr(icol)
               enddo

               ! copy out FAR taumol InOuts
               if (do_FAR) then
                  idx = gicol_cld(cols:cole)
                  gtaumol_age(idx) = tmage     (1:ncol)
                  gsflxzen (:,idx) = sflxzen (:,1:ncol)
                  gssi     (:,idx) = ssi     (:,1:ncol)
                  gtaur  (:,:,idx) = taur  (:,:,1:ncol)
                  gtaug  (:,:,idx) = taug  (:,:,1:ncol)
                 Rgcldymc(:,:,idx) = merge(1.,0.,cldymc(:,:,1:ncol))
                  gtaucmc(:,:,idx) = taucmc(:,:,1:ncol)
                  gssacmc(:,:,idx) = ssacmc(:,:,1:ncol)
                  gasmcmc(:,:,idx) = asmcmc(:,:,1:ncol)
                  gtaormc(:,:,idx) = taormc(:,:,1:ncol)
               end if

            endif  ! clear/cloudy

            deallocate(idx,__STAT__)
            call MAPL_TimerOff(MAPL,"RRTMG_PART",__RC__)

         enddo  ! over partitions

      enddo  ! outer loop (cc) over clear then cloudy columns

      ! If the user requests 'normalized' fluxes, divide
      ! the fluxes by the solar constant times coszen
      ! MAT This requires only lit points passed in

      if (normFlx == 1) then

         gswdflx_at_top(:) = max(gswdflx(:,nlay+1),1e-7)

         do ilev = 1,nlay+1
            gswuflxc(:,ilev) = gswuflxc(:,ilev) / gswdflx_at_top(:)
            gswdflxc(:,ilev) = gswdflxc(:,ilev) / gswdflx_at_top(:)
            gswuflx (:,ilev) = gswuflx (:,ilev) / gswdflx_at_top(:)
            gswdflx (:,ilev) = gswdflx (:,ilev) / gswdflx_at_top(:)
         enddo

         gnirr(:) = gnirr(:) / gswdflx_at_top(:)
         gnirf(:) = gnirf(:) / gswdflx_at_top(:)
         gparr(:) = gparr(:) / gswdflx_at_top(:)
         gparf(:) = gparf(:) / gswdflx_at_top(:)
         guvrr(:) = guvrr(:) / gswdflx_at_top(:)
         guvrf(:) = guvrf(:) / gswdflx_at_top(:)

      endif

      _RETURN(_SUCCESS)
   end subroutine rrtmg_sw_sub


   !-----------------------------------------------------------------------
   real function earth_sun(idn)
   !-----------------------------------------------------------------------
   !
   !  Purpose: Function to calculate the correction factor of Earth's orbit
   !  for current day of the year
   ! 
   !  idn        : Day of the year
   !  earth_sun  : square of the ratio of mean to actual Earth-Sun distance
   !-----------------------------------------------------------------------

      use rrsw_con, only : pi

      integer, intent(in) :: idn

      real :: gamma

      gamma = 2. * pi * (idn-1)/365. 

      ! Use Iqbal's equation 1.2.1

      earth_sun = 1.000110  + .034221  * cos(gamma) + .001289  * sin(gamma) + &
                   .000719  * cos(2. *gamma) + .000077  * sin(2. *gamma)

   end function earth_sun

end module rrtmg_sw_rad


