#include "MAPL_Generic.h"
   
module rad_utils

  use ESMF
  use MAPL

  implicit none
  private

  public :: Tbr_from_band_flux

contains

  ! estimate brightness temperature from a band flux
  subroutine Tbr_from_band_flux(IM, JM, Fband_, wn1, wn2, Tbr_, RC)

    ! input arguments
    integer, intent(in ) :: IM, JM
    real,    intent(in ) :: Fband_(IM,JM) ! band flux [W/m2]
    real,    intent(in ) :: wn1, wn2      ! bounds of band [m-1]
 
    ! output arguments
    real,    intent(out) :: Tbr_(IM,JM)   ! brightness temp [K]
 
    ! error code
    integer, optional, intent(out) :: RC
 
    ! fundamental constants
    double precision, parameter :: h  = 6.626070040d-34  ! Plancks constant         [J.s]
    double precision, parameter :: c  = 2.99792458d8     ! Speed of light in vacuum [m/s]
    double precision, parameter :: kB = 1.38064852d-23   ! Boltzmann constant       [J/K]
    double precision, parameter :: pi = MAPL_PI_R8
 
    ! other constants
    double precision, parameter :: alT = h * c / kB
    double precision, parameter :: bigS = 2.0d0 * kB**4 * pi / (h**3 * c**2)
    double precision, parameter :: bigC = 2.0d0 * h * c**2
 
    ! locals
    integer :: STATUS
    double precision, dimension(IM,JM) :: Fband, Tbr, Bmean
    real :: wnMid
 
    if (present(RC)) RC = ESMF_SUCCESS
    _ASSERT(wn2 > wn1,'band wavenumber bounds mis-ordered!')
 
    ! calculations done in double precision
    Fband = dble(Fband_)
 
    ! first guess Tbr from narrow band approximation ...
    ! (1) estimate mean Planck function for a narrow band
    Bmean = Fband / (pi * (wn2 - wn1))
    ! (2) invert Planck function for temp at mid-point wavenumber
    wnMid = (wn1 + wn2) / 2.0d0
    call invert_Planck_for_T(IM, JM, Bmean, wnMid, bigC, alT, Tbr, __RC__)
 
    ! now refine with a wide band esimate
    ! PMN: Iterative routine not ready for prime time
    !      Produces erroneously large Tbr in cloudy regions
    !call Tbr_wide_band(IM, JM, Fband, wn1, wn2, bigS, alT, Tbr, __RC__)
 
    ! put output back in real
    where (Tbr > 0.0d0)
      Tbr_ = real(Tbr)
    elsewhere
      Tbr_ = MAPL_UNDEF
    endwhere

  end subroutine Tbr_from_band_flux

  ! invert Planck function for temperature
  subroutine invert_Planck_for_T(IM, JM, Bwn, wn, bigC, alT, T, RC)

    ! input arguments
    integer,          intent(in ) :: IM, JM
    double precision, intent(in ) :: Bwn(IM,JM)  ! PlanckFn(wavenumber)
    real,             intent(in ) :: wn          ! wavenumber [m-1]
    double precision, intent(in ) :: bigC, alT   ! necessary constants

    ! output arguments
    double precision, intent(out) :: T(IM,JM)    ! temperature [K]

    ! error code
    integer, optional, intent(out) :: RC

    ! error checking
    if (present(RC)) RC = ESMF_SUCCESS
    _ASSERT(wn > 0.,'non-positive wavenumber!')

    ! invert Planck function for temp
    where (Bwn > 0.0d0)
      T = (alT * wn) / log((bigC * wn**3) / Bwn + 1.0d0)
    elsewhere
      T = 0.0d0
    endwhere

  end subroutine invert_Planck_for_T

  ! Tbr from wide band approximation
  subroutine Tbr_wide_band(IM, JM, Fband, wn1, wn2, bigS, alT, Tbr, RC)

    ! input arguments
    integer,          intent(in   ) :: IM, JM
    double precision, intent(in   ) :: Fband(IM,JM)  ! band flux [W/m2]
    real,             intent(in   ) :: wn1, wn2      ! bounds of band [m-1]
    double precision, intent(in   ) :: bigS, alT     ! necessary constant
 
    ! Tbr inputs first guess and outputs better estimate
    double precision, intent(inout) :: Tbr(IM,JM)    ! brightness temp [K]
 
    ! error code
    integer, optional, intent(out) :: RC
 
    ! number of iterations for wide band estimate (converges slowly)
    integer, parameter :: Nits = 16
 
    ! locals
    integer :: n
    real    :: alTwn1, alTwn2
 
    ! error checking
    if (present(RC)) RC = ESMF_SUCCESS
    _ASSERT(Nits >= 1,'must have at least one iteration!')
 
    ! iterate from first guess Tbr to better estimate
    alTwn1 = alT * wn1
    alTwn2 = alT * wn2
    do n = 1, Nits
      where (Tbr > 0.0d0) &
        Tbr = ( Fband / (bigS * (Tfunc(alTwn1/Tbr) - Tfunc(alTwn2/Tbr))) ) ** 0.25d0
    end do

  end subroutine Tbr_wide_band

  elemental double precision function Tfunc(x)

    double precision, intent(in) :: x

    ! maximum number of terms in series (converges quickly)
    integer, parameter :: nmax = 4

    ! locals
    integer :: n, n2, n3
    double precision :: emx, cx0, cx1, cx2, cx3, zn

    ! setup
    emx = exp(-x)
    cx0 = 6.0d0
    cx1 = 6.0d0 * x
    cx2 = 3.0d0 * x**2
    cx3 =         x**3

    ! do at least 1st order
    Tfunc = (cx3 + cx2 + cx1 + cx0) * emx
    if (nmax <= 1) return

    ! higher orders
    zn = emx
    do n = 2, nmax
      n2 = n * n
      n3 = n * n2
      zn = zn * emx
      Tfunc = Tfunc + (cx3 + cx2/n + cx1/n2 + cx0/n3) * zn / n
    end do

  end function Tfunc

end module rad_utils
