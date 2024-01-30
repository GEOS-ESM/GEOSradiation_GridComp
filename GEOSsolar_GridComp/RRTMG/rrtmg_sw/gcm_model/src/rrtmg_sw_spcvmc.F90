!pmn: table lookup creates cache issues cf cost of fn evaluation
!pmn: and must have at some point been removed ... consider
!pmn: removing LW table lookup as well.

#include "MAPL_Generic.h"

module rrtmg_sw_spcvmc

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2009, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

   ! ------- Modules -------

   use ESMF
   use MAPL

   use parrrsw, only : nbndsw, ngptsw, jpband
   use rrsw_tbl, only : od_lo
   use rrsw_wvn, only : ngb
   use rrtmg_sw_taumol, only: taumol_sw

   implicit none

contains

   ! ---------------------------------------------------------------------------
   subroutine spcvmc_sw (MAPL, &
      cc, pncol, ncol, nlay, &
      palbd, palbp, &
      pcldymc, ptaucmc, pasycmc, pomgcmc, ptaormc, &
      pltaormc, pltaucmc, pitaormc, pitaucmc, &
      ptaua, pasya, pomga, prmu0, adjflux, &
      isolvar, svar_f, svar_s, svar_i, &
      svar_f_bnd, svar_s_bnd, svar_i_bnd, &
      laytrop, jp, jt, jt1, &
      colch4, colco2, colh2o, colmol, colo2, colo3, &
      fac00, fac01, fac10, fac11, &
      cloudLM, cloudMH, & 
      selffac, selffrac, indself, forfac, forfrac, indfor, &
      pbbfd, pbbfu, pbbcd, pbbcu, puvfd, puvcd, pnifd, pnicd, &
      pbbfddir, pbbcddir, puvfddir, puvcddir, pnifddir, pnicddir, &
      znirr, znirf, zparr, zparf, zuvrr, zuvrf, fndsbnd, &
      zcotdtp, zcotdhp, zcotdmp, zcotdlp, &
      zcotntp, zcotnhp, zcotnmp, zcotnlp, &
      zcdsdtp, zcdsdhp, zcdsdmp, zcdsdlp, &
      zcdsntp, zcdsnhp, zcdsnmp, zcdsnlp, &
      zcotldtp, zcotldhp, zcotldmp, zcotldlp, &
      zcotlntp, zcotlnhp, zcotlnmp, zcotlnlp, &
      zcdsldtp, zcdsldhp, zcdsldmp, zcdsldlp, &
      zcdslntp, zcdslnhp, zcdslnmp, zcdslnlp, &
      zcotidtp, zcotidhp, zcotidmp, zcotidlp, &
      zcotintp, zcotinhp, zcotinmp, zcotinlp, &
      zcdsidtp, zcdsidhp, zcdsidmp, zcdsidlp, &
      zcdsintp, zcdsinhp, zcdsinmp, zcdsinlp, &
      do_drfband, zdrband, zdfband, &
      RC)
   ! ---------------------------------------------------------------------------
   !
   ! Purpose: Contains spectral loop to compute the shortwave radiative fluxes, 
   !          using the two-stream method of H. Barker and McICA, the Monte-Carlo
   !          Independent Column Approximation, for the representation of 
   !          sub-grid cloud variability (i.e. cloud overlap).
   !
   ! Interface:  *spcvmc_sw* is called from *rrtmg_sw.F90* or rrtmg_sw.1col.F90*
   !
   ! Method:
   !    Adapted from two-stream model of H. Barker;
   !    Two-stream model options (selected with kmodts in rrtmg_sw_reftra.F90):
   !        1: Eddington, 2: PIFM, Zdunkowski et al., 3: discret ordinates
   !
   ! Modifications:
   !
   ! Original: H. Barker
   ! Revision: Merge with RRTMG_SW: J.-J.Morcrette, ECMWF, Feb 2003
   ! Revision: Add adjustment for Earth/Sun distance : MJIacono, AER, Oct 2003
   ! Revision: Bug fix for use of PALBP and PALBD: MJIacono, AER, Nov 2003
   ! Revision: Bug fix to apply delta scaling to clear sky: AER, Dec 2004
   ! Revision: Code modified so that delta scaling is not done in cloudy profiles
   !           if routine cldprop is used; delta scaling can be applied by swithcing
   !           code below if cldprop is not used to get cloud properties. 
   !           AER, Jan 2005
   ! Revision: Modified to use McICA: MJIacono, AER, Nov 2005
   ! Revision: Uniform formatting for RRTMG: MJIacono, AER, Jul 2006 
   ! Revision: Use exponential lookup table for transmittance: MJIacono, AER, 
   !           Aug 2007 
   ! Revision: Added ztau[lmht]p: PMNorris, GMAO, at some point
   ! Revision: Added zdrband, zdfband for OBIO support: PMNorris, GMAO, Nov 2022
   ! Revision: separate phase tracking for taormc, taucmc: PMNorris, Jan 2024
   !
   ! ------------------------------------------------------------------

      ! ------- Input -------

      type(MAPL_MetaComp), pointer, intent(inout) :: MAPL

      integer, intent(in) :: pncol, ncol, cc
      integer, intent(in) :: nlay
      integer, intent(in) :: laytrop (pncol)

      integer, intent(in) :: jp  (nlay,pncol) 
      integer, intent(in) :: jt  (nlay,pncol) 
      integer, intent(in) :: jt1 (nlay,pncol) 

      real, intent(in) :: adjflux(jpband)             ! Earth/Sun distance adjustment

      ! Solar variability
      integer, intent(in) :: isolvar                  ! Flag for solar variability method
      real, intent(in) :: svar_f                      ! Solar variability facular multiplier
      real, intent(in) :: svar_s                      ! Solar variability sunspot multiplier
      real, intent(in) :: svar_i                      ! Solar variability baseline irradiance multiplier
      real, intent(in) :: svar_f_bnd (jpband)         ! Solar variability facular multiplier (by band)
      real, intent(in) :: svar_s_bnd (jpband)         ! Solar variability sunspot multiplier (by band)
      real, intent(in) :: svar_i_bnd (jpband)         ! Solar variability baseline irradiance multiplier (by band)

      real, intent(in) :: palbd (nbndsw,pncol)        ! surface albedo (diffuse)
      real, intent(in) :: palbp (nbndsw,pncol)        ! surface albedo (direct)
      real, intent(in) :: prmu0        (pncol)        ! cosine of solar zenith angle

      ! McICA cloud presence and optical properties:
      ! These are only used, and therefore only need to be defined, for cloudy gridcolumns
      ! (and hence potentially cloudy subcolumns) (cc==2) not for clear gridcolumns (cc==1)
      logical, intent(in) :: pcldymc (nlay,ngptsw,pncol)  ! cloudy or not? [mcica]
      real,    intent(in) :: ptaucmc (nlay,ngptsw,pncol)  ! cloud optical depth [mcica]
      real,    intent(in) :: pasycmc (nlay,ngptsw,pncol)  ! cloud asymmetry parameter [mcica]
      real,    intent(in) :: pomgcmc (nlay,ngptsw,pncol)  ! cloud single scattering albedo [mcica]
      real,    intent(in) :: ptaormc (nlay,ngptsw,pncol)  ! cloud optical depth, non-delta scaled [mcica]
   
      real,    intent(in) :: pltaormc(nlay,ngptsw,pncol)  ! liq cloud optical depth, non-delta scaled [mcica]
      real,    intent(in) :: pltaucmc(nlay,ngptsw,pncol)  ! liq cloud optical depth [mcica]
      real,    intent(in) :: pitaormc(nlay,ngptsw,pncol)  ! ice cloud optical depth, non-delta scaled [mcica]
      real,    intent(in) :: pitaucmc(nlay,ngptsw,pncol)  ! ice cloud optical depth [mcica]

      real, intent(in) :: ptaua (nlay,nbndsw,pncol)  ! aerosol optical depth
      real, intent(in) :: pasya (nlay,nbndsw,pncol)  ! aerosol asymmetry parameter
      real, intent(in) :: pomga (nlay,nbndsw,pncol)  ! aerosol single scattering albedo
                                                               
      real, intent(in) :: colh2o (nlay,pncol) 
      real, intent(in) :: colco2 (nlay,pncol) 
      real, intent(in) :: colch4 (nlay,pncol) 
      real, intent(in) :: colo3  (nlay,pncol) 
      real, intent(in) :: colo2  (nlay,pncol) 
      real, intent(in) :: colmol (nlay,pncol) 

      ! continuum interpolation coefficients
      integer, intent(in) :: indself  (nlay,pncol)
      integer, intent(in) :: indfor   (nlay,pncol)
      real,    intent(in) :: selffac  (nlay,pncol)
      real,    intent(in) :: selffrac (nlay,pncol)
      real,    intent(in) :: forfac   (nlay,pncol)
      real,    intent(in) :: forfrac  (nlay,pncol)

      ! pressure and temperature interpolation coefficients
      real, intent(in), dimension (nlay,pncol) &
         :: fac00, fac01, fac10, fac11

      ! pressure super-layer interface levels for optical thicknesses
      integer, intent(in) :: cloudLM  ! Low-mid
      integer, intent(in) :: cloudMH  ! Mid-high

      logical, intent(in) :: do_drfband  ! Compute zdrband, zdfband?

      ! ------- Output -------

      real, intent(out) :: pbbcd    (nlay+1,pncol) 
      real, intent(out) :: pbbcu    (nlay+1,pncol) 
      real, intent(out) :: pbbfd    (nlay+1,pncol) 
      real, intent(out) :: pbbfu    (nlay+1,pncol) 
      real, intent(out) :: pbbfddir (nlay+1,pncol) 
      real, intent(out) :: pbbcddir (nlay+1,pncol) 

      real, intent(out) :: puvcd    (nlay+1,pncol) 
      real, intent(out) :: puvfd    (nlay+1,pncol) 
      real, intent(out) :: puvcddir (nlay+1,pncol) 
      real, intent(out) :: puvfddir (nlay+1,pncol) 

      real, intent(out) :: pnicd    (nlay+1,pncol) 
      real, intent(out) :: pnifd    (nlay+1,pncol) 
      real, intent(out) :: pnicddir (nlay+1,pncol) 
      real, intent(out) :: pnifddir (nlay+1,pncol) 
      
      ! surface band fluxes (direct and TOTAL)
      real, intent(out), dimension(pncol) :: &
         znirr, znirf, zparr, zparf, zuvrr, zuvrf

      ! net downwelling flux @ sfc in bands (all-sky and diffuse+direct)
      real, intent(out) :: fndsbnd (pncol,nbndsw)

      ! in-cloud PAR optical thicknesses
      real, intent(out), dimension (pncol) :: &
         zcotdtp, zcotdhp, zcotdmp, zcotdlp, &
         zcotntp, zcotnhp, zcotnmp, zcotnlp
      real, intent(out), dimension (pncol) :: &
         zcdsdtp, zcdsdhp, zcdsdmp, zcdsdlp, &
         zcdsntp, zcdsnhp, zcdsnmp, zcdsnlp
      real, intent(out), dimension (pncol) :: &
         zcotldtp, zcotldhp, zcotldmp, zcotldlp, &
         zcotlntp, zcotlnhp, zcotlnmp, zcotlnlp
      real, intent(out), dimension (pncol) :: &
         zcdsldtp, zcdsldhp, zcdsldmp, zcdsldlp, &
         zcdslntp, zcdslnhp, zcdslnmp, zcdslnlp
      real, intent(out), dimension (pncol) :: &
         zcotidtp, zcotidhp, zcotidmp, zcotidlp, &
         zcotintp, zcotinhp, zcotinmp, zcotinlp
      real, intent(out), dimension (pncol) :: &
         zcdsidtp, zcdsidhp, zcdsidmp, zcdsidlp, &
         zcdsintp, zcdsinhp, zcdsinmp, zcdsinlp

      ! Surface downwelling direct and diffuse fluxes (W/m2)
      !    in each band (all-sky): Only filled if (do_drfband).
      real, intent(out), dimension (pncol,nbndsw) :: zdrband, zdfband

      integer, intent(out), optional :: RC  ! return code

      ! ------- Local -------

      integer :: icol
      integer :: jk, ikl
      integer :: iw, jb, ibm

      real :: zf, zwf, zincflx, wgt
      real :: staotp, staohp, staomp, staolp
      real :: stautp, stauhp, staump, staulp
      real :: sltaotp, sltaohp, sltaomp, sltaolp
      real :: sltautp, sltauhp, sltaump, sltaulp
      real :: sitaotp, sitaohp, sitaomp, sitaolp
      real :: sitautp, sitauhp, sitaump, sitaulp

      real :: zgco   (nlay,ngptsw,pncol)
      real :: zomco  (nlay,ngptsw,pncol)  
      real :: ztauo  (nlay,ngptsw,pncol)  

      real :: zdbt   (nlay,  ngptsw,pncol)
      real :: ztdbt  (nlay+1,ngptsw,pncol)   

      real :: zfd    (nlay+1,ngptsw,pncol)
      real :: zfu    (nlay+1,ngptsw,pncol)   

      real :: zref   (nlay+1,ngptsw,pncol)  ! direct beam reflectivity
      real :: zrefd  (nlay+1,ngptsw,pncol)  ! diffuse     reflectivity
      real :: ztra   (nlay+1,ngptsw,pncol)  ! direct beam transmissivity
      real :: ztrad  (nlay+1,ngptsw,pncol)  ! diffuse     transmissivity

      real :: ztaur  (nlay,ngptsw,pncol)
      real :: ztaug  (nlay,ngptsw,pncol) 

      real :: zsflxzen (ngptsw,pncol)
      real :: ssi      (ngptsw,pncol)

      integer :: STATUS  ! for MAPL error reporting

      ! ------------------------------------------------------------------

      ! zero output accumulators
      pbbcd    = 0. 
      pbbcu    = 0. 
      pbbfd    = 0. 
      pbbfu    = 0. 
      pbbcddir = 0. 
      pbbfddir = 0. 
      puvcd    = 0. 
      puvfd    = 0. 
      puvcddir = 0. 
      puvfddir = 0. 
      pnicd    = 0. 
      pnifd    = 0. 
      pnicddir = 0. 
      pnifddir = 0.
      znirr    = 0.
      znirf    = 0.
      zparr    = 0.
      zparf    = 0.
      zuvrr    = 0.
      zuvrf    = 0.
      fndsbnd  = 0.
      if (do_drfband) then
         zdrband  = 0.
         zdfband  = 0.
      end if

      ! Calculate the optical depths for gaseous absorption and Rayleigh scattering     
      call MAPL_TimerOn(MAPL,"---RRTMG_TAUMOL",__RC__)
      call taumol_sw( &
         pncol, ncol, nlay, &
         colh2o, colco2, colch4, colo2, colo3, colmol, &
         laytrop, jp, jt, jt1, fac00, fac01, fac10, fac11, &
         selffac, selffrac, indself, forfac, forfrac, indfor, &
         isolvar, svar_f, svar_s, svar_i, svar_f_bnd, svar_s_bnd, svar_i_bnd, &
         ssi, zsflxzen, ztaug, ztaur)
      call MAPL_TimerOff(MAPL,"---RRTMG_TAUMOL",__RC__)

      ! Set fixed boundary values.
      ! The sfc (jk=nlay+1) zref[d] & ztra[d] never change from these.
      ! The TOA (jk=1) ztdbt likewise never change.
      do icol = 1,ncol
         do iw = 1,ngptsw
            jb = ngb(iw)  ! SW band: jb = 16:29
            ibm = jb-15   !   => ibm = 1:14

            ! Surface reflectivities / transmissivities
            zref  (nlay+1,iw,icol) = palbp(ibm,icol) 
            zrefd (nlay+1,iw,icol) = palbd(ibm,icol) 
            ztra  (nlay+1,iw,icol) = 0. 
            ztrad (nlay+1,iw,icol) = 0. 
           
            ! TOA direct beam    
            ztdbt (1,iw,icol) = 1. 
    
         end do
      end do

      ! Delta-scaled clear-sky optical properties
      do icol = 1,ncol
         do iw = 1,ngptsw
            jb = ngb(iw)
            ibm = jb-15
            do jk = 1,nlay
               ikl = nlay+1-jk

               ! Clear-sky optical parameters including aerosols
               ! Gas: ssa=0, g=0, Rayleigh: ssa=1, g=0
               ztauo(jk,iw,icol) = ztaur(ikl,iw,icol) + ztaug(ikl,iw,icol) + ptaua(ikl,ibm,icol)      
               zomco(jk,iw,icol) = ztaur(ikl,iw,icol) + ptaua(ikl,ibm,icol) * pomga(ikl,ibm,icol)
               zgco (jk,iw,icol) = (pasya(ikl,ibm,icol) * pomga(ikl,ibm,icol) * ptaua(ikl,ibm,icol)) &
                                   / zomco(jk,iw,icol)   
               zomco(jk,iw,icol) = zomco(jk,iw,icol) / ztauo(jk,iw,icol)
               
               ! Delta-scaling with f = g**2
               zf = zgco(jk,iw,icol) ** 2
               zwf = zomco(jk,iw,icol) * zf
               ztauo(jk,iw,icol) = (1. - zwf) * ztauo(jk,iw,icol)
               zomco(jk,iw,icol) = (zomco(jk,iw,icol) - zwf) / (1. - zwf)
               zgco (jk,iw,icol) = (zgco (jk,iw,icol) - zf ) / (1. - zf )

            enddo    
         end do
      end do

      ! Clear-sky reflectivities / transmissivities
      ! note: pcldymc may not be defined here but the
      !       last arg .false. means it is not used anyway.
      call MAPL_TimerOn(MAPL,"---RRTMG_REFTRA",__RC__)
      call reftra_sw (pncol, ncol, nlay, &
                      pcldymc, zgco, prmu0, ztauo, zomco, &
                      zref, zrefd, ztra, ztrad, .false.)
      call MAPL_TimerOff(MAPL,"---RRTMG_REFTRA",__RC__)

      ! Clear-sky direct beam transmittance        
      do icol = 1,ncol
         do iw = 1,ngptsw
            do jk = 1,nlay
               zdbt(jk,iw,icol) = exp(-ztauo(jk,iw,icol) / prmu0(icol))
               ztdbt(jk+1,iw,icol) = zdbt(jk,iw,icol) * ztdbt(jk,iw,icol)  
            enddo          
         end do
      end do

      ! Vertical quadrature for clear-sky fluxes
      call MAPL_TimerOn(MAPL,"---RRTMG_VRTQDR",__RC__)
      call vrtqdr_sw(pncol, ncol, nlay, &
                     zref, zrefd, ztra, ztrad, &
                     zdbt, ztdbt, &
                     zfd, zfu)
      call MAPL_TimerOff(MAPL,"---RRTMG_VRTQDR",__RC__)

      ! Band integration for clear cases      
      do icol = 1,ncol
         do iw = 1,ngptsw
            jb = ngb(iw)
            ibm = jb-15

            ! Adjust incoming flux for Earth/Sun distance and zenith angle
            if (isolvar < 0) then
               ! No solar variability and no solar cycle
               zincflx = adjflux(jb) * zsflxzen(iw,icol) * prmu0(icol)           
            else
               ! Solar variability with averaged or specified solar cycle
               zincflx = adjflux(jb) * ssi     (iw,icol) * prmu0(icol)           
            endif

            do ikl = 1,nlay+1
               jk = nlay+2-ikl

               ! Accumulate spectral fluxes over whole spectrum  
               pbbcu   (ikl,icol) = pbbcu   (ikl,icol) + zincflx * zfu  (jk,iw,icol)  
               pbbcd   (ikl,icol) = pbbcd   (ikl,icol) + zincflx * zfd  (jk,iw,icol)  
               pbbcddir(ikl,icol) = pbbcddir(ikl,icol) + zincflx * ztdbt(jk,iw,icol)  
              
               ! Band fluxes
               if (ibm >= 10 .and. ibm <= 13) then
                  ! UV/Visible
                  puvcd   (ikl,icol)  = puvcd   (ikl,icol)  + zincflx * zfd  (jk,iw,icol)  
                  puvcddir(ikl,icol)  = puvcddir(ikl,icol)  + zincflx * ztdbt(jk,iw,icol)  
               else if (ibm == 14 .or. ibm <= 9) then  
                  ! Near-IR
                  pnicd   (ikl,icol) = pnicd   (ikl,icol) + zincflx * zfd  (jk,iw,icol)  
                  pnicddir(ikl,icol) = pnicddir(ikl,icol) + zincflx * ztdbt(jk,iw,icol)  
               endif 

            enddo  ! layer loop
         enddo  ! spectral loop
      enddo  ! column loop

      !!!!!!!!!!!!!!!!
      !! END CLEAR  !!
      !!!!!!!!!!!!!!!!

      if (cc == 2) then

         ! Add in cloud optical properties (which are already delta-scaled)
         ! (uses the pre-calculated clear-sky delta-scaled zomco and zgco)
         do icol = 1,ncol
            do iw = 1,ngptsw
               jb = ngb(iw)
               ibm = jb-15
               do jk = 1,nlay
                  ikl = nlay+1-jk

                  ! only cloudy gridboxes need to add in cloud opt props
                  if (pcldymc(ikl,iw,icol)) then
                     zgco (jk,iw,icol) = &
                        ztauo  (jk, iw,icol) * zomco  (jk, iw,icol) * zgco   (jk, iw,icol) + &
                        ptaucmc(ikl,iw,icol) * pomgcmc(ikl,iw,icol) * pasycmc(ikl,iw,icol)
                     zomco(jk,iw,icol) = &
                        ztauo  (jk,iw, icol) * zomco  (jk, iw,icol) + &
                        ptaucmc(ikl,iw,icol) * pomgcmc(ikl,iw,icol)
                     ztauo(jk,iw,icol) = &
                        ztauo  (jk, iw,icol) + &
                        ptaucmc(ikl,iw,icol) 
                     zgco (jk,iw,icol) = zgco (jk,iw,icol) / zomco(jk,iw,icol)
                     zomco(jk,iw,icol) = zomco(jk,iw,icol) / ztauo(jk,iw,icol)
                  end if
               
               enddo    
            end do
         end do

         ! Update reflectivities / transmissivities for cloudy cells only
         ! note: since cc==2 here pcldymc is defined
         call MAPL_TimerOn(MAPL,"---RRTMG_REFTRA",__RC__)
         call reftra_sw (pncol, ncol, nlay, &
                         pcldymc, zgco, prmu0, ztauo, zomco, &
                         zref, zrefd, ztra, ztrad, .true.)
         call MAPL_TimerOff(MAPL,"---RRTMG_REFTRA",__RC__)

         ! Recalculate direct transmission
         do icol = 1,ncol
            do iw = 1,ngptsw
               do jk = 1,nlay
                  ikl = nlay+1-jk 

                  ! cloudy ztauo has been updated so recalculate
                  if (pcldymc(ikl,iw,icol)) &
                     zdbt(jk,iw,icol) = exp(-ztauo(jk,iw,icol) / prmu0(icol))            
                  ztdbt(jk+1,iw,icol) = zdbt(jk,iw,icol) * ztdbt(jk,iw,icol)

               enddo          
            end do
         end do

         ! Vertical quadrature for total-sky fluxes
         call MAPL_TimerOn(MAPL,"---RRTMG_VRTQDR",__RC__)
         call vrtqdr_sw(pncol, ncol, nlay, &
                        zref, zrefd, ztra, ztrad, &
                        zdbt, ztdbt, &
                        zfd, zfu)
         call MAPL_TimerOff(MAPL,"---RRTMG_VRTQDR",__RC__)

         ! Upwelling and downwelling fluxes at levels
         !   Two-stream calculations go from top to bottom; 
         !   layer indexing is reversed to go bottom to top for output arrays

         do icol = 1,ncol
            do iw = 1,ngptsw
               jb = ngb(iw)
               ibm = jb-15

               ! Adjust incoming flux for Earth/Sun distance and zenith angle
               if (isolvar < 0) then
                  ! No solar variability and no solar cycle
                  zincflx = adjflux(jb) * zsflxzen(iw,icol) * prmu0(icol)           
               else
                  ! Solar variability with averaged or specified solar cycle
                  zincflx = adjflux(jb) * ssi     (iw,icol) * prmu0(icol)           
               endif

               do ikl = 1,nlay+1
                  jk = nlay+2-ikl

                  ! Accumulate spectral fluxes over whole spectrum  
                  pbbfu   (ikl,icol)  = pbbfu   (ikl,icol) + zincflx * zfu  (jk,iw,icol)  
                  pbbfd   (ikl,icol)  = pbbfd   (ikl,icol) + zincflx * zfd  (jk,iw,icol)              
                  pbbfddir(ikl,icol)  = pbbfddir(ikl,icol) + zincflx * ztdbt(jk,iw,icol)  

                  ! Band fluxes
                  if (ibm >= 10 .and. ibm <= 13) then
                     ! UV/Visible
                     puvfd   (ikl,icol) = puvfd   (ikl,icol) + zincflx * zfd  (jk,iw,icol)  
                     puvfddir(ikl,icol) = puvfddir(ikl,icol) + zincflx * ztdbt(jk,iw,icol)  
                  else if (ibm == 14 .or. ibm <= 9) then  
                     ! Near-IR
                     pnifd   (ikl,icol) = pnifd   (ikl,icol) + zincflx * zfd  (jk,iw,icol)  
                     pnifddir(ikl,icol) = pnifddir(ikl,icol) + zincflx * ztdbt(jk,iw,icol)  
                  endif

               enddo  ! layer loop
            enddo  ! spectral loop
         enddo  ! column loop

      else

         ! cc == 1 so clear so total-sky == clear-sky
         pbbfu    = pbbcu
         pbbfd    = pbbcd
         pbbfddir = pbbcddir
         puvfd    = puvcd
         puvfddir = puvcddir
         pnifd    = pnicd
         pnifddir = pnicddir

      end if

      ! surface band fluxes
      do icol = 1,ncol
         do iw = 1,ngptsw
            jb = ngb(iw)
            ibm = jb - 15

            ! Adjust incoming flux for Earth/Sun distance and zenith angle
            if (isolvar < 0) then
               ! No solar variability and no solar cycle
               zincflx = adjflux(jb) * zsflxzen(iw,icol) * prmu0(icol)           
            else
               ! Solar variability with averaged or specified solar cycle
               zincflx = adjflux(jb) * ssi     (iw,icol) * prmu0(icol)           
            endif

            ! Band fluxes
            if (ibm == 14 .or. ibm <= 8) then
               ! near-IR
               znirr(icol) = znirr(icol) + zincflx * ztdbt(nlay+1,iw,icol)  ! Direct flux
               znirf(icol) = znirf(icol) + zincflx * zfd  (nlay+1,iw,icol)  ! Total flux
            else if (ibm >= 10 .and. ibm <= 11) then
               ! Photosynthetically active radiation (PAR)
               zparr(icol) = zparr(icol) + zincflx * ztdbt(nlay+1,iw,icol)  ! Direct flux
               zparf(icol) = zparf(icol) + zincflx * zfd  (nlay+1,iw,icol)  ! Total flux
            else if (ibm >= 12 .and. ibm <= 13) then
               ! UV
               zuvrr(icol) = zuvrr(icol) + zincflx * ztdbt(nlay+1,iw,icol)  ! Direct flux
               zuvrf(icol) = zuvrf(icol) + zincflx * zfd  (nlay+1,iw,icol)  ! Total flux
            else if ( ibm==9) then
               ! Band 9: half PAR, half NIR
               zparr(icol) = zparr(icol) + 0.5 * zincflx * ztdbt(nlay+1,iw,icol)  ! Direct flux
               zparf(icol) = zparf(icol) + 0.5 * zincflx * zfd  (nlay+1,iw,icol)  ! Total flux
               znirr(icol) = znirr(icol) + 0.5 * zincflx * ztdbt(nlay+1,iw,icol)  ! Direct flux
               znirf(icol) = znirf(icol) + 0.5 * zincflx * zfd  (nlay+1,iw,icol)  ! Total flux
            endif

            fndsbnd(icol,ibm) = fndsbnd(icol,ibm) + &
               zincflx * (zfd (nlay+1,iw,icol) - zfu (nlay+1,iw,icol))

            if (do_drfband) then
               zdrband(icol,ibm) = zdrband(icol,ibm) + zincflx * ztdbt(nlay+1,iw,icol)  ! direct
               zdfband(icol,ibm) = zdfband(icol,ibm) + zincflx * zfd  (nlay+1,iw,icol)  ! total
            end if

         end do
      enddo                    

      ! convert from total to diffuse only
      if (do_drfband) zdfband = zdfband - zdrband

      ! diagnostic in-cloud optical thicknesses in PAR super-band
      ! (weighted across and within bands by TOA incident flux)
      ! -------------------------------------------------------
      zcotdtp = 0.; zcotntp = 0.
      zcotdhp = 0.; zcotnhp = 0.
      zcotdmp = 0.; zcotnmp = 0.
      zcotdlp = 0.; zcotnlp = 0.
      zcdsdtp = 0.; zcdsntp = 0.
      zcdsdhp = 0.; zcdsnhp = 0.
      zcdsdmp = 0.; zcdsnmp = 0.
      zcdsdlp = 0.; zcdsnlp = 0.
      zcotldtp = 0.; zcotlntp = 0.
      zcotldhp = 0.; zcotlnhp = 0.
      zcotldmp = 0.; zcotlnmp = 0.
      zcotldlp = 0.; zcotlnlp = 0.
      zcdsldtp = 0.; zcdslntp = 0.
      zcdsldhp = 0.; zcdslnhp = 0.
      zcdsldmp = 0.; zcdslnmp = 0.
      zcdsldlp = 0.; zcdslnlp = 0.
      zcotidtp = 0.; zcotintp = 0.
      zcotidhp = 0.; zcotinhp = 0.
      zcotidmp = 0.; zcotinmp = 0.
      zcotidlp = 0.; zcotinlp = 0.
      zcdsidtp = 0.; zcdsintp = 0.
      zcdsidhp = 0.; zcdsinhp = 0.
      zcdsidmp = 0.; zcdsinmp = 0.
      zcdsidlp = 0.; zcdsinlp = 0.

      ! can only be non-zero for potentially cloudy columns
      if (cc == 2) then

         do icol = 1,ncol

            do iw = 1,ngptsw
               jb = ngb(iw)
               ibm = jb - 15
   
               ! band weights
               if (ibm >= 10 .and. ibm <= 11) then
                  ! Photosynthetically active radiation (PAR)
                  wgt = 1.0
               else if (ibm == 9) then
                  ! half PAR, half NIR
                  wgt = 0.5
               else
                  ! no contribution to PAR
                  cycle
               end if
   
               ! TOA flux weighting (adjustment for zenith angle
               ! not needed since normalized for each icol)
               if (isolvar < 0) then
                  zincflx = adjflux(jb) * zsflxzen(iw,icol)
               else
                  zincflx = adjflux(jb) * ssi     (iw,icol)
               endif
               wgt = wgt * zincflx

               ! low pressure layer
               staolp = sum(ptaormc(1:cloudLM,iw,icol),dim=1)
               if (staolp > 0.) then
                  zcotnlp(icol) = zcotnlp(icol) + wgt * staolp
                  zcotdlp(icol) = zcotdlp(icol) + wgt
               end if
               staulp = sum(ptaucmc(1:cloudLM,iw,icol),dim=1)
               if (staulp > 0.) then
                  zcdsnlp(icol) = zcdsnlp(icol) + wgt * staulp
                  zcdsdlp(icol) = zcdsdlp(icol) + wgt
               end if
               sltaolp = sum(pltaormc(1:cloudLM,iw,icol),dim=1)
               if (sltaolp > 0.) then
                  zcotlnlp(icol) = zcotlnlp(icol) + wgt * sltaolp
                  zcotldlp(icol) = zcotldlp(icol) + wgt
               end if
               sltaulp = sum(pltaucmc(1:cloudLM,iw,icol),dim=1)
               if (sltaulp > 0.) then
                  zcdslnlp(icol) = zcdslnlp(icol) + wgt * sltaulp
                  zcdsldlp(icol) = zcdsldlp(icol) + wgt
               end if
               sitaolp = sum(pitaormc(1:cloudLM,iw,icol),dim=1)
               if (sitaolp > 0.) then
                  zcotinlp(icol) = zcotinlp(icol) + wgt * sitaolp
                  zcotidlp(icol) = zcotidlp(icol) + wgt
               end if
               sitaulp = sum(pitaucmc(1:cloudLM,iw,icol),dim=1)
               if (sitaulp > 0.) then
                  zcdsinlp(icol) = zcdsinlp(icol) + wgt * sitaulp
                  zcdsidlp(icol) = zcdsidlp(icol) + wgt
               end if

               ! mid pressure layer
               staomp = sum(ptaormc(cloudLM+1:cloudMH,iw,icol),dim=1)
               if (staomp > 0.) then
                  zcotnmp(icol) = zcotnmp(icol) + wgt * staomp
                  zcotdmp(icol) = zcotdmp(icol) + wgt
               end if
               staump = sum(ptaucmc(cloudLM+1:cloudMH,iw,icol),dim=1)
               if (staump > 0.) then
                  zcdsnmp(icol) = zcdsnmp(icol) + wgt * staump
                  zcdsdmp(icol) = zcdsdmp(icol) + wgt
               end if
               sltaomp = sum(pltaormc(cloudLM+1:cloudMH,iw,icol),dim=1)
               if (sltaomp > 0.) then
                  zcotlnmp(icol) = zcotlnmp(icol) + wgt * sltaomp
                  zcotldmp(icol) = zcotldmp(icol) + wgt
               end if
               sltaump = sum(pltaucmc(cloudLM+1:cloudMH,iw,icol),dim=1)
               if (sltaump > 0.) then
                  zcdslnmp(icol) = zcdslnmp(icol) + wgt * sltaump
                  zcdsldmp(icol) = zcdsldmp(icol) + wgt
               end if
               sitaomp = sum(pitaormc(cloudLM+1:cloudMH,iw,icol),dim=1)
               if (sitaomp > 0.) then
                  zcotinmp(icol) = zcotinmp(icol) + wgt * sitaomp
                  zcotidmp(icol) = zcotidmp(icol) + wgt
               end if
               sitaump = sum(pitaucmc(cloudLM+1:cloudMH,iw,icol),dim=1)
               if (sitaump > 0.) then
                  zcdsinmp(icol) = zcdsinmp(icol) + wgt * sitaump
                  zcdsidmp(icol) = zcdsidmp(icol) + wgt
               end if

               ! high pressure layer
               staohp = sum(ptaormc(cloudMH+1:nlay,iw,icol),dim=1)
               if (staohp > 0.) then
                  zcotnhp(icol) = zcotnhp(icol) + wgt * staohp
                  zcotdhp(icol) = zcotdhp(icol) + wgt
               end if
               stauhp = sum(ptaucmc(cloudMH+1:nlay,iw,icol),dim=1)
               if (stauhp > 0.) then
                  zcdsnhp(icol) = zcdsnhp(icol) + wgt * stauhp
                  zcdsdhp(icol) = zcdsdhp(icol) + wgt
               end if
               sltaohp = sum(pltaormc(cloudMH+1:nlay,iw,icol),dim=1)
               if (sltaohp > 0.) then
                  zcotlnhp(icol) = zcotlnhp(icol) + wgt * sltaohp
                  zcotldhp(icol) = zcotldhp(icol) + wgt
               end if
               sltauhp = sum(pltaucmc(cloudMH+1:nlay,iw,icol),dim=1)
               if (sltauhp > 0.) then
                  zcdslnhp(icol) = zcdslnhp(icol) + wgt * sltauhp
                  zcdsldhp(icol) = zcdsldhp(icol) + wgt
               end if
               sitaohp = sum(pitaormc(cloudMH+1:nlay,iw,icol),dim=1)
               if (sitaohp > 0.) then
                  zcotinhp(icol) = zcotinhp(icol) + wgt * sitaohp
                  zcotidhp(icol) = zcotidhp(icol) + wgt
               end if
               sitauhp = sum(pitaucmc(cloudMH+1:nlay,iw,icol),dim=1)
               if (sitauhp > 0.) then
                  zcdsinhp(icol) = zcdsinhp(icol) + wgt * sitauhp
                  zcdsidhp(icol) = zcdsidhp(icol) + wgt
               end if

               ! whole subcolumn
               staotp = staolp + staomp + staohp
               if (staotp > 0.) then
                  zcotntp(icol) = zcotntp(icol) + wgt * staotp
                  zcotdtp(icol) = zcotdtp(icol) + wgt
               end if
               stautp = staulp + staump + stauhp
               if (stautp > 0.) then
                  zcdsntp(icol) = zcdsntp(icol) + wgt * stautp
                  zcdsdtp(icol) = zcdsdtp(icol) + wgt
               end if
               sltaotp = sltaolp + sltaomp + sltaohp
               if (sltaotp > 0.) then
                  zcotlntp(icol) = zcotlntp(icol) + wgt * sltaotp
                  zcotldtp(icol) = zcotldtp(icol) + wgt
               end if
               sltautp = sltaulp + sltaump + sltauhp
               if (sltautp > 0.) then
                  zcdslntp(icol) = zcdslntp(icol) + wgt * sltautp
                  zcdsldtp(icol) = zcdsldtp(icol) + wgt
               end if
               sitaotp = sitaolp + sitaomp + sitaohp
               if (sitaotp > 0.) then
                  zcotintp(icol) = zcotintp(icol) + wgt * sitaotp
                  zcotidtp(icol) = zcotidtp(icol) + wgt
               end if
               sitautp = sitaulp + sitaump + sitauhp
               if (sitautp > 0.) then
                  zcdsintp(icol) = zcdsintp(icol) + wgt * sitautp
                  zcdsidtp(icol) = zcdsidtp(icol) + wgt
               end if

            end do ! iw
         end do  ! icol
      end if  ! cc==2

      _RETURN(_SUCCESS)
   end subroutine spcvmc_sw

   ! --------------------------------------------------------------------
   subroutine reftra_sw(pncol, ncol, nlay, &
                        cloudy, pgg, prmuzl, ptau, pw, &
                        pref, prefd, ptra, ptrad, &
                        update_cloudy_cells_only)
   ! --------------------------------------------------------------------
   ! Purpose: computes the reflectivity and transmissivity of a clear or 
   !   cloudy layers using a choice of various approximations.
   !
   ! inputs:
   !      cloudy  = binary cloud presence flag
   !      pgg     = assymetry factor
   !      prmuzl  = cosine solar zenith angle
   !      ptau    = optical thickness
   !      pw      = single scattering albedo
   !
   ! outputs:
   !      pref    : collimated beam reflectivity
   !      prefd   : diffuse beam reflectivity 
   !      ptra    : collimated beam transmissivity
   !      ptrad   : diffuse beam transmissivity
   !
   ! notes:
   !   1. This routine is intended to be run twice. First with clear-only
   !      optical parameters and update_cloudy_cells_only false. This gives
   !      the clear-sky reflectivities/transmissivities for every cell.
   !      The second time through, it is run with cloud-adjusted optical
   !      parameters and update_cloudy_cells_only true, which only updates
   !      cells that actually have clouds. The clear-cells retain their
   !      former clear outputs via the intent(inout) "outputs".
   !   2. The cloudy flag is only used if update_cloudy_cells_only true,
   !      so only necessary to define it in that case.
   !
   ! Method:
   !      standard delta-eddington, p.i.f.m., or d.o.m. layer calculations.
   !      kmodts  = 1 eddington (joseph et al., 1976)
   !              = 2 pifm (zdunkowski et al., 1980)	<--
   !              = 3 discrete ordinates (liou, 1973)
   !
   ! Modifications:
   ! Original: J-JMorcrette, ECMWF, Feb 2003
   ! Revised for F90 reformatting: MJIacono, AER, Jul 2006
   ! Revised to add exponential lookup table: MJIacono, AER, Aug 2007
   ! Cleaned up and modified as per above note: PMNorris, GMAO, Aug 2021
   !
   ! ------------------------------------------------------------------

      ! ------- Input -------

      integer, intent (in) :: pncol  ! dimensioned num of gridcols
      integer, intent (in) :: ncol   ! actual number of gridcols
      integer, intent (in) :: nlay   ! number of model layers


      logical, intent(in) :: cloudy  (nlay,ngptsw,pncol)  ! cloudy or not?

      real,    intent(in) :: pgg     (nlay,ngptsw,pncol)  ! asymmetry parameter
      real,    intent(in) :: ptau    (nlay,ngptsw,pncol)  ! optical depth
      real,    intent(in) :: pw      (nlay,ngptsw,pncol)  ! single scattering albedo 
      real,    intent(in) :: prmuzl              (pncol)  ! cosine of solar zenith angle

      logical, intent(in) :: update_cloudy_cells_only

      ! ------- Output -------

      real, intent(inout) :: pref  (nlay+1,ngptsw,pncol)  ! direct beam reflectivity
      real, intent(inout) :: prefd (nlay+1,ngptsw,pncol)  ! diffuse     reflectivity
      real, intent(inout) :: ptra  (nlay+1,ngptsw,pncol)  ! direct beam transmissivity
      real, intent(inout) :: ptrad (nlay+1,ngptsw,pncol)  ! diffuse     transmissivity

      ! ------- Local -------

      integer :: jk, kmodts
      integer :: icol, iw
      logical :: update

      real :: za, za1, za2
      real :: zbeta, zdend, zdenr, zdent
      real :: ze1, ze2, zem1, zem2, zemm, zep1, zep2
      real :: zg, zg3, zgamma1, zgamma2, zgamma3, zgamma4, zgt
      real :: zr1, zr2, zr3, zr4, zr5
      real :: zrk, zrk2, zrkg, zrm1, zrp, zrp1, zrpp
      real :: zsr3, zt1, zt2, zt3, zt4, zt5, zto1
      real :: zw, zwcrit, zwo, prmuz

      real, parameter :: eps = 1.e-08 

      ! MAT These are 8-byte versions of zw, zg, and zwo. This is done
      ! MAT to avoid a divide-by-zero in the zwo calculation below. More
      ! MAT information below.
      ! MAT NOTE: This is not an official fix, just a patch to allow work
      ! MAT       for now.

      real*8 :: zw8, zg8, zwo8

      ! ------------------------------------------------------------------

      zsr3 = sqrt(3.)
      zwcrit = 0.9999995 
      kmodts = 2
      
      do icol = 1,ncol
         prmuz = prmuzl(icol)
         do iw = 1,ngptsw
            do jk = 1,nlay

               if (update_cloudy_cells_only) then
                  update = cloudy(nlay+1-jk,iw,icol)
               else
                  update = .true.
               end if
               if (update) then

                  zto1 = ptau(jk,iw,icol)  
                  zw   = pw  (jk,iw,icol)  
                  zg   = pgg (jk,iw,icol)    

                  ! MAT Move zw and zg into 8-byte reals to avoid
                  ! MAT divide-by-zero in zwo calculation below

                  zw8 = zw
                  zg8 = zg

                  ! General two-stream expressions

                  zg3 = 3. * zg
           
                  zgamma1 = (8. - zw * (5. + zg3)) * 0.25 
                  zgamma2 = 3. * (zw * (1. - zg)) * 0.25 
                  zgamma3 = (2. - zg3 * prmuz) * 0.25 
                  zgamma4 = 1. - zgamma3
    
                  ! Recompute original s.s.a. to test for conservative solution

                  ! MAT The issue with this is as follows. A column occurs for
                  ! MAT which zw = 0.0249478 and zg = 0.503158. If you
                  ! MAT calculate the denominator of zwo, you get a value of
                  ! MAT -0.000000064995 which, apparently, is then flushed to
                  ! MAT zero by the compiler. Then zwo = zw / 0. and bam.
                  ! MAT To avoid this loss of precision, we calculate zwo
                  ! MAT in 8-byte reals.
                  ! MAT
                  ! MAT Original code
                  ! MAT zwo= zw / (1.  - (1.  - zw) * (zg / (1.  - zg))**2)
                  ! MAT
                  ! MAT New code
                  zwo8 = zw8 / (1.d0 - (1.d0 - zw8) * (zg8 / (1.d0 - zg8))**2)

                  ! MAT Put zwo8 into a 4-byte real
                  zwo = zwo8

                  ! END MAT EDITS
    
                  if (zwo >= zwcrit) then
                     ! Conservative scattering

                     za  = zgamma1 * prmuz 
                     za1 = za - zgamma3
                     zgt = zgamma1 * zto1
        
                     ! Homogeneous reflectance and transmittance,
                     ! collimated beam

                     ze1 = min(zto1/prmuz, 500.)
                     ze2 = exp(-ze1)
                     pref(jk,iw,icol) = (zgt - za1 * (1. - ze2)) / (1. + zgt)
                     ptra(jk,iw,icol) = 1. - pref(jk,iw,icol)  

                     ! isotropic incidence

                     prefd(jk,iw,icol) = zgt / (1. + zgt)
                     ptrad(jk,iw,icol) = 1. - prefd(jk,iw,icol)          

                     ! This is applied for consistency between total (delta-scaled) and direct (unscaled) 
                     ! calculations at very low optical depths (tau < 1.e-4) when the exponential lookup
                     ! table returns a transmittance of 1.
                     if (ze2 == 1.) then 
                        pref (jk,iw,icol) = 0. 
                        ptra (jk,iw,icol) = 1. 
                        prefd(jk,iw,icol) = 0. 
                        ptrad(jk,iw,icol) = 1. 
                     endif

                  else
                     ! Non-conservative scattering

                     za1 = zgamma1 * zgamma4 + zgamma2 * zgamma3
                     za2 = zgamma1 * zgamma3 + zgamma2 * zgamma4
                     zrk = sqrt(zgamma1**2 - zgamma2**2)
                     zrp = zrk * prmuz               
                     zrp1 = 1. + zrp
                     zrm1 = 1. - zrp
                     zrk2 = 2. * zrk
                     zrpp = 1. - zrp*zrp
                     zrkg = zrk + zgamma1
                     zr1  = zrm1 * (za2 + zrk * zgamma3)
                     zr2  = zrp1 * (za2 - zrk * zgamma3)
                     zr3  = zrk2 * (zgamma3 - za2 * prmuz)
                     zr4  = zrpp * zrkg
                     zr5  = zrpp * (zrk - zgamma1)
                     zt1  = zrp1 * (za1 + zrk * zgamma4)
                     zt2  = zrm1 * (za1 - zrk * zgamma4)
                     zt3  = zrk2 * (zgamma4 + za1 * prmuz)
                     zt4  = zr4
                     zt5  = zr5

                     ! mji - reformulated code to avoid potential floating point exceptions
                     !               zbeta = - zr5 / zr4
                     zbeta = (zgamma1 - zrk) / zrkg
        
                     ! Homogeneous reflectance and transmittance

                     ze1 = min(zrk * zto1, 5.)
                     ze2 = min(zto1 / prmuz, 5.)

                     ! exponential transmittance (or expansion of exp()for low tau)
                     ! pmn: run with od_lo approx commented out was a little slower,
                     !      so have kept it, even though pure exp() is more accurate
                     !      and avoids a discontinuity.
                     if (ze1 <= od_lo) then 
                        zem1 = 1. - ze1 + 0.5 * ze1 * ze1
                     else
                        zem1 = exp(-ze1)
                     endif
                     zep1 = 1. / zem1
                     if (ze2 <= od_lo) then 
                        zem2 = 1. - ze2 + 0.5 * ze2 * ze2
                     else
                        zem2 = exp(-ze2)
                     endif
                     zep2 = 1. / zem2

                     zdenr = zr4*zep1 + zr5*zem1
                     zdent = zt4*zep1 + zt5*zem1
                     if (zdenr >= -eps .and. zdenr <= eps) then
                        pref(jk,iw,icol) = eps
                        ptra(jk,iw,icol) = zem2
                     else 
                        pref(jk,iw,icol) = zw * (zr1*zep1 - zr2*zem1 - zr3*zem2) / zdenr
                        ptra(jk,iw,icol) = zem2 - zem2 * zw * (zt1*zep1 - zt2*zem1 - zt3*zep2) / zdent
                     endif

                     ! diffuse beam

                     zemm = zem1*zem1
                     zdend = 1. / ((1. - zbeta * zemm) * zrkg)
                     prefd(jk,iw,icol) = zgamma2 * (1. - zemm) * zdend
                     ptrad(jk,iw,icol) = zrk2 * zem1 * zdend

                  endif
               endif         

            end do  
         end do
      end do

   end subroutine reftra_sw
                           

   ! --------------------------------------------------------------------------
   subroutine vrtqdr_sw(pncol, ncol, nlay, &
                        pref, prefd, ptra, ptrad, &
                        pdbt, ptdbt, &
                        pfd, pfu)
   ! --------------------------------------------------------------------------
   ! Purpose: This routine performs the vertical quadrature integration
   !
   ! Interface:  *vrtqdr_sw* is called from *spcvrt_sw* and *spcvmc_sw*
   !
   ! Modifications.
   ! 
   ! Original: H. Barker
   ! Revision: Integrated with rrtmg_sw, J.-J. Morcrette, ECMWF, Oct 2002
   ! Revision: Reformatted for consistency with rrtmg_lw: MJIacono, AER, Jul 2006
   ! Revision: Tidy up, reordering for speed, etc.: PM Norris, GMAO, Aug 2021
   ! Revision: Add code explanations: PM Norris, GMAO, Dec 2022.
   !
   !-----------------------------------------------------------------------

      ! ----- Input -----

      integer, intent (in) :: pncol                   ! dimensioned num of gridcols
      integer, intent (in) :: ncol                    ! actual number of gridcols
      integer, intent (in) :: nlay                    ! number of model layers
    
      ! Layer reflectivities and transmittivities (1..nlay). The nlay+1 values
      ! of pref[d] have albedo values discussed below. Note: the direct/diffuse
      ! nomenclature refers to the SOURCE prior to reflection/transmission. Of
      ! course a diffuse source remains diffuse in its reflection/tranmission,
      ! but a direct beam source will, in general, produce both direct beam and
      ! diffuse components in transmission (& a diffuse-only reflected component).
      real, intent(in) :: pref (nlay+1,ngptsw,pncol)  ! direct beam reflectivity
      real, intent(in) :: prefd(nlay+1,ngptsw,pncol)  ! diffuse reflectivity
      real, intent(in) :: ptra (nlay+1,ngptsw,pncol)  ! direct beam transmissivity
      real, intent(in) :: ptrad(nlay+1,ngptsw,pncol)  ! diffuse transmissivity

      ! Note: these describe the transmission of the direct beam as a direct beam.
      real, intent(in) :: pdbt (nlay,  ngptsw,pncol)  ! lyr mean dir beam transmittance
      real, intent(in) :: ptdbt(nlay+1,ngptsw,pncol)  ! total direct beam transmittance

      ! ----- Output -----
      ! unadjusted for earth/sun distance or zenith angle

      real, intent(out) :: pfd(nlay+1,ngptsw,pncol)   ! downwelling flux (W/m2)
      real, intent(out) :: pfu(nlay+1,ngptsw,pncol)   ! upwelling   flux (W/m2)
    
      ! ----- Local -----

      integer :: ikp, ikx, jk, icol, iw
      real :: zreflect, zreflectj

      real :: ztdn  (nlay+1,ngptsw,pncol)
      real :: prup  (nlay+1,ngptsw,pncol)  
      real :: prupd (nlay+1,ngptsw,pncol)  
      real :: prdnd (nlay+1,ngptsw,pncol)  
     
      do icol = 1,ncol
         do iw = 1,ngptsw
      
            ! The nlay+1 level of prup/prupd require palbp/palbd which
            ! are already available in fixed nlay+1 level of pref/prefd.
            ! Explanation: prup[d](nlay+1) are reflected fractions of
            ! downelling direct[diffuse] radiation due to the surface alone.

            prup (nlay+1,iw,icol) = pref (nlay+1,iw,icol)
            prupd(nlay+1,iw,icol) = prefd(nlay+1,iw,icol)

            ! Link lowest layer with surface.
            ! Explanation: prup[d](nlay) are reflected fractions of downwelling
            ! direct[diffuse] radiation due to the lowest model layer and the
            ! surface combined. The first term in each, pref[d](nlay), accounts
            ! for the reflection from the lowest model layer alone, before
            ! reaching the surface. Any additioinal reflected radiation must
            ! involve a surface reflection (which results in diffuse reflected
            ! radiation regardless of type of incoming radiation) and needs a
            ! a final transmission through the lowest model layer to make it
            ! out. This is the ptrad(nlay) factor common to both prup[d](nlay).
            ! The factor ptrad(nlay) * prefd(nlay+1) in prupd(nlay) is the 
            ! necessary transmission of the diffuse incoming radiation to the
            ! surface and its reflection by the surface. The corresponding
            ! factor for prup(nlay) is more complicated:
            ! (ptra(nlay) - pdbt(nlay)) * prefd(nlay+1) + pdbt(nlay) * pref(nlay+1)
            ! Here pdbt(nlay) is the fraction of the incoming direct beam trans-
            ! mitted to the surface as a direct beam, while ptra(nlay) - pdbt(nlay)
            ! is the fraction of the incoming direct beam transmitted to the surface
            ! as diffuse radiation. Each is multiplied by the appropriate surface
            ! reflection, pref(nlay+1) and prefd(nlay+1), respectively. Finally,
            ! AFTER reflection by the surface, and BEFORE final transmission thru
            ! the model layer, the reflected diffuse radiation can undergo ANY 
            ! NUMBER of RE-reflections by the layer, prefd(nlay), and subsequent
            ! reflections, prefd(nlay+1), by the surface. This gives rise to a
            ! geometric series that converges to the zreflect term.

            zreflect = 1. / (1. - prefd(nlay+1,iw,icol) * prefd(nlay,iw,icol))
            prup(nlay,iw,icol) = pref(nlay,iw,icol) + (ptrad(nlay,iw,icol) * &
               ((ptra(nlay,iw,icol) - pdbt(nlay,iw,icol)) * prefd(nlay+1,iw,icol) + &
               pdbt(nlay,iw,icol) * pref(nlay+1,iw,icol))) * zreflect
            prupd(nlay,iw,icol) = prefd(nlay,iw,icol) + ptrad(nlay,iw,icol) * &
               ptrad(nlay,iw,icol) * prefd(nlay+1,iw,icol) * zreflect

          end do
      end do
      
      ! Pass from bottom to top 
      ! Explantion: prup[d](ikx) are reflected fractions of downwelling direct[diffuse]
      ! radiation due to layer ikx, all those below it, and the surface, combined.
      ! This is achieved by isolating the effects of layer ikx alone, and the cmbined
      ! effects of the lower layers (ikp .. nlay) and the surface. This proceeds by 
      ! iteration, adding one layer ikx at a time and depending on the previous 
      ! prup[d](ikp). Otherwise it is completely analagopus to the explation above,
      ! but with prup[d](ikp) for the layers/surface below ikx repacing prefd(nlay+1)
      ! for the surface in the above formulas/explanation.
      ! Note that the the top of leyer ikx is level ikx, which is another way of
      ! thinking about the index ikx in these definitions and the formula below.
      ! PMN: possible to combine this loop with one above?

      do icol = 1,ncol
         do iw = 1,ngptsw
            do jk = 1,nlay-1

               ikp = nlay+1-jk                       
               ikx = ikp-1
               zreflectj = 1. / (1. - prupd(ikp,iw,icol) * prefd(ikx,iw,icol))
               prup(ikx,iw,icol) = pref(ikx,iw,icol) + (ptrad(ikx,iw,icol) * &
                  ((ptra(ikx,iw,icol) - pdbt(ikx,iw,icol)) * prupd(ikp,iw,icol) + &
                  pdbt(ikx,iw,icol) * prup(ikp,iw,icol))) * zreflectj
               prupd(ikx,iw,icol) = prefd(ikx,iw,icol) + ptrad(ikx,iw,icol) * &
                  ptrad(ikx,iw,icol) * prupd(ikp,iw,icol) * zreflectj

            enddo
         end do
      end do

      ! Now work from top to bottom calculating:
      ! ztdn(ikp), the transmision fraction of the incoming top of
      !    model direct beam flux through the bottom of layer ikp-1,
      !    as either direct beam or diffuse radiation.
      ! prdnd(ikp), the combined reflection fraction of upwelling
      !    diffuse radiation downward by layers ikp-1 and above.
      ! Note that the the bottom of leyer ikp-1 is level ikp, which
      ! is another way of thinking about the index ikp in these
      ! definitions and the formula below.

      do icol = 1,ncol
         do iw = 1,ngptsw

            ! Upper boundary conditions

            ztdn (1,iw,icol) = 1. 
            prdnd(1,iw,icol) = 0. 
            ztdn (2,iw,icol) = ptra (1,iw,icol)  
            prdnd(2,iw,icol) = prefd(1,iw,icol)  

         end do
      end do
      
      do icol = 1,ncol
         do iw = 1,ngptsw

            ! Pass from top to bottom
            ! Explanation: prdnd is simple in terms of the definition above and
            ! the concepts discussed in previous explanations. For ztdn, first
            ! note the use of ptdbt(jk) (cf. pdbt), which is the total direct-
            ! beam-remaining transmission from the top of the model thru LEVEL jk.
            ! This is one of only two routes --- not involving reflection(s) from
            ! layer jk below --- for radiation to arrive at level jk. The other
            ! is as diffuse radiation ztdn(jk) - ptdbt(jk).
            !   In the first case, the direct-beam source, ptdbt(jk), is either
            ! directly transmitted by layer jk, the term ptdbt(jk) * ptra(jk),
            ! or else it is indirectly transmitted as diffuse, factor ptrad(jk),
            ! after:
            !     (1) a double-reflection pref(jk) * prdnd(jk), namely a direct-
            !   source reflection by layer jk, pref(jk), converting it to diffuse,
            !   and a re-reflection down from layers above jk, prdnd(jk);
            !     (2) followed by zero or more all-diffuse double-reflections
            !   prefd(jk) * prdnd(jk), as encapsulated in zreflect.
            ! Note that it is because of the differences between the reflection
            ! and transmission of a layer for direct-source vs. diffuse-source,
            ! namely the difference between pref and prefd and ptra and ptrad,
            ! that this first case appears as it does and cannot be simplified.
            !   In the second case, the diffuse fraction ztdn(jk) - ptdbt(jk)
            ! is transmitted after zero or more double-reflections, hence the
            ! final fraction, (ztdn(jk) - ptdbt(jk)) * zreflect * ptrad(jk).
            ! Note that the radiation is always diffuse for this case, so the
            ! complexity of the first case is not required.

            do jk = 2,nlay
               ikp = jk+1

               zreflect = 1. / (1. - prefd(jk,iw,icol) * prdnd(jk,iw,icol))
               ztdn(ikp,iw,icol) = ptdbt(jk,iw,icol) * ptra(jk,iw,icol) + &
                  (ptrad(jk,iw,icol) * ((ztdn(jk,iw,icol) - ptdbt(jk,iw,icol)) + &
                  ptdbt(jk,iw,icol) * pref(jk,iw,icol) * prdnd(jk,iw,icol))) * zreflect
               prdnd(ikp,iw,icol) = prefd(jk,iw,icol) + ptrad(jk,iw,icol) * &
                  ptrad(jk,iw,icol) * prdnd(jk,iw,icol) * zreflect

            enddo
         end do
      end do
    
      ! Up and down-welling fluxes at levels

      do icol = 1,ncol
         do iw = 1,ngptsw
            do jk = 1,nlay+1
               zreflect = 1. / (1. - prdnd(jk,iw,icol) * prupd(jk,iw,icol))
               pfu(jk,iw,icol) = (ptdbt(jk,iw,icol) * prup(jk,iw,icol) + &
                  (ztdn(jk,iw,icol) - ptdbt(jk,iw,icol)) * prupd(jk,iw,icol)) * zreflect
               pfd(jk,iw,icol) = ptdbt(jk,iw,icol) + (ztdn(jk,iw,icol) - ptdbt(jk,iw,icol) + &
                  ptdbt(jk,iw,icol) * prup(jk,iw,icol) * prdnd(jk,iw,icol)) * zreflect
            enddo
         end do
      end do
      
   end subroutine vrtqdr_sw

end module rrtmg_sw_spcvmc
