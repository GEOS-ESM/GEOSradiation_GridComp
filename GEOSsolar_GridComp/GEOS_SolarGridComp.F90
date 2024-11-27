#include "MAPL_Generic.h"
#define LIN2_ARG1(VAR,I,J,FINT) (VAR(I,J) + FINT * (VAR(I+1,J)-VAR(I,J)))

! ==============================================================================
! Note: the SOLAR_RADVAL compile time flag (enabled with the ENABLE_SOLAR_RADVAL 
! CMake option) is used to select solar diagnostic features which are generally
! more advanced than what a regular user will need and mainly for use by the
! the radiation code development team. They are chosen by compile time flag 
! because they bloat the restart state and may also incur other computational
! costs that are not warranted under normal (non-development) use.
! ==============================================================================

module GEOS_SolarGridCompMod

!=============================================================================
!BOP

! !MODULE: GEOS_SolarGridCompMod -- Computes solar radiation fluxes in a cloudy atmosphere

! /*
! !DESCRIPTION:
!
! {\tt GEOS\_SolarGridCompMod} is an ESMF/MAPL gridded component that performs
! a broadband calculation of shortwave radiative fluxes for use as a solar
! radiation parameterization in atmospheric models on a sphere. \newline
!
! {\em Scientific Basis:} The radiative transfer calculation is based on the
! M-D Chou shortwave parameterization. The basic reference for the scheme is:
! Chou and Suarez 1999: A Solar Radiation Parameterization for Atmospheric
! Studies, NASA-TM-1999-104606, Vol 15. An updated version of this report
! can be found in SolarDoc.pdf in this directory. \newline
!
! The parameterization treats direct and diffuse fluxes of solar radiation
! in eight spectral bands:
! \begin{verbatim}
!        in the uv region :
!           index  1 for the 0.225-0.285 micron band
!           index  2 for the 0.175-0.225;0.285-0.300 micron band
!           index  3 for the 0.300-0.325 micron band
!           index  4 for the 0.325-0.4 micron band
!        in the par region :
!           index  5 for the 0.4-0.690 micron band
!        in the infrared region :
!           index  6 for the 0.690-1.220 micron band
!           index  7 for the 1.220-2.270 micron band
!           index  8 for the 2.270-3.850 micron band
! \end{verbatim}
! It includes gaseous absorption due to water vapor, ozone, CO$_2$, and
! molecular oxygen and the effects of molecular scattering, as well as
! multiple scattering due to clouds and aerosols. \newline
!
! It allows clouds to occur in any layer and horizontal cloud cover fractions
! must be specified for all layers; clear layers simply have a zero fraction.
! Vertically, the layers are assumed to be filled by cloud. To simplify the
! treatment of cloud effects, the model layers, are grouped into three super
! layers. Effective cloud properties are then parameterized by assuming that
! clouds are maximally overlapped within the super layers and randomly over-
! lapped between the super layers. The optical properties of cloud particles
! depend on the liquid, ice, and rain mixing ratios, as well as on spatially
! dependent effective radii for the three species. These are all inputs to
! the Gridded Component. \newline
!
! The parameterization can include the effects of an arbitrary number of
! aerosol species. Aerosol optical thickness, single-scattering albedo, and
! asymmetry factor must be determined as functions of height and spectral
! band for each species. \newline
!
! {\em Code Implementation:} \newline
!
! {\tt GEOS\_SolarGridCompMod} is an encapsulation of Chou's plug-compatible
! SORAD Fortran routine in a MAPL/ESMF gridded component (GC). It follows the
! standard rules for an ESMF/MAPL GCs. It operates on the ESMF grid that
! appears in the gridded component. This grid must be present in the GC and
! properly initialized before Initialize is called. The only restrictions on
! the grid are that it be 3-dimensional with two horizontal and one vertical
! dimension and with only the horizontal dimensions decomposed. The vertical
! dimension is also assumed to the the third dimension of the Fortran arrays
! and is indexed from the top down. No particular vertical coordinate is
! assumed, rather the 3-dimensional field of air pressure at the layer
! interfaces is a required Import. \newline
!
! This module contains only SetServices and Run methods. The Initialize
! and Finalize methods being defaulted to the MAPL\_Generic versions. The
! SetServices method is the only public entity. There are no public types
! or data. \newline
!
! The contents of the Import, Export, and Internal States are explicitly
! described in SetServices and in tables in this documentation. All quantities
! in these states are in either ESMF Fields or Bundles, and all share a common
! grid---the ESMF grid in the gridded component at the time Initialize (in this
! case, MAPL\_GenericInitialize) was called. All outputs appearing in the Export
! state are optional and are filled only if they have been allocated. All filled
! Exports are valid for the time interval on the GC's clock when the run method
! is invoked. Imports can be from either an instantaneous or a time-averaged
! state of the atmosphere. All Imports are read-only; none are Friendly. Most
! imports are simple ESMF Fields containing 2- or 3-dimensional quantities,
! such as temperature and humidity, needed in the flux calculation. Non-cloud
! aerosol amounts are the exception; they appear in an ESMF Bundle. \newline
!
! The net (+ve downward) fluxes on the Export state are defined at the layer
! interfaces, which are indexed from the top of the atmosphere (L=0) to the
! surface. Incident fluxes at the surface also appear in the Export state;
! these are separated into direct (beam) and diffuse fluxes for three spectral
! bands (uv, par, nir), as defined in the table above. \newline
!
! The full transfer calculation is done infrequently and its results kept in
! the Internal state. The frequency of full calculations is controlled by an
! alarm whose interval can be set from a value in the configuration and whose
! origin is taken as the beginning of the run. For the full calculations, solar
! fluxes are computed based on mean zenith angles averaged over sun positions
! for a given period (the long interval, which can be specified in the config-
! uration) beyond the current time on the input clock. On every call to the Run
! method, whatever the state of the alarm that controls the full calculation,
! the sun's position is updated to the mean position for the clock's current
! interval and fluxes are updated based on normalized fluxes computed during
! the previous full transfer calculation, but using the TOA insolation for the
! current time on the clock. Because of this intermittent scheme, checkpoint-
! restart sequences are seamless only when interrupted at the time of the full
! calculation. \newline
!
! The calculation relies in MAPL's Astronomy layer, which in turn assumes that
! the ESMF grid can be queried for latitude and longitude coordinates. \newline
!
! {\em Configuration:} \newline
!
! Like all MAPL GCs, {\tt GEOS\_SolarGridCompMod} assumes that the configuration
! in the ESMF GC is open and treats it as an environment from which it can {\em
! at any time} read control information. It uses MAPL rules for scanning this
! configuration.
!\begin{verbatim}
!
! VARIABLE             DESCRIPTION           UNITS      DEFAULT   NOTES
!
! RUN_DT:              Short time interval   (seconds)  none
! DT:                  Long time interval    (seconds)  RUN_DT
! AVGR:                Averaging interval    (seconds)  DT
! PRS_LOW_MID_CLOUDS:  Interface pressure    (Pa)       70000.
!                      between the low and
!                      middle cloud layers
! PRS_MID_HIGH_CLOUDS: Interface pressure    (Pa)       40000.
!                      between the high and
!                      middle cloud layers
! SOLAR_CONSTANT:                            (W m-2)    none      Use -1 for time-dependent values
! CO2:                 CO2 concentration     (ppmv)     none      Use -1 for time-dependent values
!
!\end{verbatim}
!
! !BUGS:
!
!\begin{itemize}
!  \item Aerosol properties for each aerosol in the Bundle are obtained by calling a global
!  method (Get\_AeroOptProp) that must recognize the aerosol by its Field name in the Bundle.
!  This is a placeholder for a scheme in which each Field carries with it a method for
!  computing its aerosol's optical properties.
!
!  \item The grid must have two horizontal dimensions and they must be the inner dimensions
!  of Fortran arrays.
!
!  \item The load-balancing relies on the grid describing a sphere. Everything works for
!  non-spherical grids but the load-balancing should be disabled and this can be done only
!  by going into the code.
!\end{itemize}
!
! */

! !USES:

  use ESMF
  use MAPL
  use gFTL_StringVector

  ! for RRTMGP
  use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp

  use soradmod, only: SORAD
  use sorad_constants, only : HK_IR_OLD, HK_UV_OLD
  use gettau, only: getvistau

  use rrtmg_sw_rad, only: rrtmg_sw
  use rrtmg_sw_init, only: rrtmg_sw_ini
  use parrrsw, only: ngptsw
  use cloud_subcol_gen, only: &
     generate_stochastic_clouds, clearCounts_threeBand

  use mo_rte_kind, only: wp

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

! !GLOBAL PARAMETERS
  INTEGER, PARAMETER :: NB_CHOU_UV  = 5 ! Number of UV bands
  INTEGER, PARAMETER :: NB_CHOU_NIR = 3 ! Number of near-IR bands
  INTEGER, PARAMETER :: NB_CHOU     = NB_CHOU_UV + NB_CHOU_NIR ! Total number of bands
  INTEGER, PARAMETER :: NB_RRTMG    = 14
  INTEGER, PARAMETER :: NB_RRTMGP   = 14
  INTEGER, PARAMETER :: NB_OBIO     = 33
  integer            :: DO_OBIO

!EOP

  ! RRTMGP internal state:
  ! Will be attached to the Gridded Component.
  ! Used to provide efficient initialization
  type ty_RRTMGP_state
    private
    logical :: initialized = .false.
    type (ty_gas_optics_rrtmgp) :: k_dist
  end type ty_RRTMGP_state

  ! Wrapper to access RRTMGP internal state
  type ty_RRTMGP_wrap
    type (ty_RRTMGP_state), pointer :: ptr => null()
  end type ty_RRTMGP_wrap

  ! ----------------------------------------------------------------
  ! For an RRTMGP forwice calculation approximating RRTMG iceflag=3:
  ! ----------------------------------------------------------------
  ! These fdlice3_rrtmgp are exactly the fdlice3 of RRTMG except they
  ! have been reordered (band 14 of RRTMG becomes band 1 of RRTMGP).
  ! The small discrepancy in the upper wavenumr edge of RRTMGP band1
  ! is ignored. Its a difference of only 80 cm-1 or 3.73 um for RRTMGP
  ! vs. 3.85 um for RRTMG). Values are unitless.
  real(wp), parameter :: fdlice3_rrtmgp(46,14) = real(reshape([&
    ! RRTMGP band 1 (820-2680 cm-1) ~ RRTMG band 14 (820-2600 cm-1)
    1.006055e-01,9.549582e-02,9.063960e-02,8.602900e-02,8.165612e-02,&
    7.751308e-02,7.359199e-02,6.988496e-02,6.638412e-02,6.308156e-02,&
    5.996942e-02,5.703979e-02,5.428481e-02,5.169657e-02,4.926719e-02,&
    4.698880e-02,4.485349e-02,4.285339e-02,4.098061e-02,3.922727e-02,&
    3.758547e-02,3.604733e-02,3.460497e-02,3.325051e-02,3.197604e-02,&
    3.077369e-02,2.963558e-02,2.855381e-02,2.752050e-02,2.652776e-02,&
    2.556772e-02,2.463247e-02,2.371415e-02,2.280485e-02,2.189670e-02,&
    2.098180e-02,2.005228e-02,1.910024e-02,1.811781e-02,1.709709e-02,&
    1.603020e-02,1.490925e-02,1.372635e-02,1.247363e-02,1.114319e-02,&
    9.727157e-03, &
    ! RRTMGP band 2 (2680-3250) ~ RRTMG band 1 (2600-3250 cm-1)
    4.959277e-02,4.685292e-02,4.426104e-02,4.181231e-02,3.950191e-02,&
    3.732500e-02,3.527675e-02,3.335235e-02,3.154697e-02,2.985578e-02,&
    2.827395e-02,2.679666e-02,2.541909e-02,2.413640e-02,2.294378e-02,&
    2.183639e-02,2.080940e-02,1.985801e-02,1.897736e-02,1.816265e-02,&
    1.740905e-02,1.671172e-02,1.606585e-02,1.546661e-02,1.490917e-02,&
    1.438870e-02,1.390038e-02,1.343939e-02,1.300089e-02,1.258006e-02,&
    1.217208e-02,1.177212e-02,1.137536e-02,1.097696e-02,1.057210e-02,&
    1.015596e-02,9.723704e-03,9.270516e-03,8.791565e-03,8.282026e-03,&
    7.737072e-03,7.151879e-03,6.521619e-03,5.841467e-03,5.106597e-03,&
    4.312183e-03, &
    ! RRTMGP band 3 = RRTMG band 2 (3250-4000 cm-1)
    5.071224e-02,5.000217e-02,4.933872e-02,4.871992e-02,4.814380e-02,&
    4.760839e-02,4.711170e-02,4.665177e-02,4.622662e-02,4.583426e-02,&
    4.547274e-02,4.514007e-02,4.483428e-02,4.455340e-02,4.429544e-02,&
    4.405844e-02,4.384041e-02,4.363939e-02,4.345340e-02,4.328047e-02,&
    4.311861e-02,4.296586e-02,4.282024e-02,4.267977e-02,4.254248e-02,&
    4.240640e-02,4.226955e-02,4.212995e-02,4.198564e-02,4.183462e-02,&
    4.167494e-02,4.150462e-02,4.132167e-02,4.112413e-02,4.091003e-02,&
    4.067737e-02,4.042420e-02,4.014854e-02,3.984840e-02,3.952183e-02,&
    3.916683e-02,3.878144e-02,3.836368e-02,3.791158e-02,3.742316e-02,&
    3.689645e-02, &
    ! RRTMGP band 4 = RRTMG band 3 (4000-4650 cm-1)
    1.062938e-01,1.065234e-01,1.067822e-01,1.070682e-01,1.073793e-01,&
    1.077137e-01,1.080693e-01,1.084442e-01,1.088364e-01,1.092439e-01,&
    1.096647e-01,1.100970e-01,1.105387e-01,1.109878e-01,1.114423e-01,&
    1.119004e-01,1.123599e-01,1.128190e-01,1.132757e-01,1.137279e-01,&
    1.141738e-01,1.146113e-01,1.150385e-01,1.154534e-01,1.158540e-01,&
    1.162383e-01,1.166045e-01,1.169504e-01,1.172741e-01,1.175738e-01,&
    1.178472e-01,1.180926e-01,1.183080e-01,1.184913e-01,1.186405e-01,&
    1.187538e-01,1.188291e-01,1.188645e-01,1.188580e-01,1.188076e-01,&
    1.187113e-01,1.185672e-01,1.183733e-01,1.181277e-01,1.178282e-01,&
    1.174731e-01, &
    ! RRTMGP band 5 = RRTMG band 4 (4650-5150 cm-1)
    1.076195e-01,1.065195e-01,1.054696e-01,1.044673e-01,1.035099e-01,&
    1.025951e-01,1.017203e-01,1.008831e-01,1.000808e-01,9.931116e-02,&
    9.857151e-02,9.785939e-02,9.717230e-02,9.650774e-02,9.586322e-02,&
    9.523623e-02,9.462427e-02,9.402484e-02,9.343544e-02,9.285358e-02,&
    9.227675e-02,9.170245e-02,9.112818e-02,9.055144e-02,8.996974e-02,&
    8.938056e-02,8.878142e-02,8.816981e-02,8.754323e-02,8.689919e-02,&
    8.623517e-02,8.554869e-02,8.483724e-02,8.409832e-02,8.332943e-02,&
    8.252807e-02,8.169175e-02,8.081795e-02,7.990419e-02,7.894796e-02,&
    7.794676e-02,7.689809e-02,7.579945e-02,7.464834e-02,7.344227e-02,&
    7.217872e-02, &
    ! RRTMGP band 6 = RRTMG band 5 (5150-6150 cm-1)
    1.119014e-01,1.122706e-01,1.126690e-01,1.130947e-01,1.135456e-01,&
    1.140199e-01,1.145154e-01,1.150302e-01,1.155623e-01,1.161096e-01,&
    1.166703e-01,1.172422e-01,1.178233e-01,1.184118e-01,1.190055e-01,&
    1.196025e-01,1.202008e-01,1.207983e-01,1.213931e-01,1.219832e-01,&
    1.225665e-01,1.231411e-01,1.237050e-01,1.242561e-01,1.247926e-01,&
    1.253122e-01,1.258132e-01,1.262934e-01,1.267509e-01,1.271836e-01,&
    1.275896e-01,1.279669e-01,1.283134e-01,1.286272e-01,1.289063e-01,&
    1.291486e-01,1.293522e-01,1.295150e-01,1.296351e-01,1.297104e-01,&
    1.297390e-01,1.297189e-01,1.296480e-01,1.295244e-01,1.293460e-01,&
    1.291109e-01, &
    ! RRTMGP band 7 = RRTMG band 6 (6150-7700 cm-1)
    1.133298e-01,1.136777e-01,1.140556e-01,1.144615e-01,1.148934e-01,&
    1.153492e-01,1.158269e-01,1.163243e-01,1.168396e-01,1.173706e-01,&
    1.179152e-01,1.184715e-01,1.190374e-01,1.196108e-01,1.201897e-01,&
    1.207720e-01,1.213558e-01,1.219389e-01,1.225194e-01,1.230951e-01,&
    1.236640e-01,1.242241e-01,1.247733e-01,1.253096e-01,1.258309e-01,&
    1.263352e-01,1.268205e-01,1.272847e-01,1.277257e-01,1.281415e-01,&
    1.285300e-01,1.288893e-01,1.292173e-01,1.295118e-01,1.297710e-01,&
    1.299927e-01,1.301748e-01,1.303154e-01,1.304124e-01,1.304637e-01,&
    1.304673e-01,1.304212e-01,1.303233e-01,1.301715e-01,1.299638e-01,&
    1.296983e-01, &
    ! RRTMGP band 8 = RRTMG band 7 (7700-8050 cm-1)
    1.145360e-01,1.153256e-01,1.161453e-01,1.169929e-01,1.178666e-01,&
    1.187641e-01,1.196835e-01,1.206227e-01,1.215796e-01,1.225522e-01,&
    1.235383e-01,1.245361e-01,1.255433e-01,1.265579e-01,1.275779e-01,&
    1.286011e-01,1.296257e-01,1.306494e-01,1.316703e-01,1.326862e-01,&
    1.336951e-01,1.346950e-01,1.356838e-01,1.366594e-01,1.376198e-01,&
    1.385629e-01,1.394866e-01,1.403889e-01,1.412678e-01,1.421212e-01,&
    1.429469e-01,1.437430e-01,1.445074e-01,1.452381e-01,1.459329e-01,&
    1.465899e-01,1.472069e-01,1.477819e-01,1.483128e-01,1.487976e-01,&
    1.492343e-01,1.496207e-01,1.499548e-01,1.502346e-01,1.504579e-01,&
    1.506227e-01, &
    ! RRTMGP band 9 = RRTMG band 8 (8050-12850 cm-1)
    1.153263e-01,1.161445e-01,1.169932e-01,1.178703e-01,1.187738e-01,&
    1.197016e-01,1.206516e-01,1.216217e-01,1.226099e-01,1.236141e-01,&
    1.246322e-01,1.256621e-01,1.267017e-01,1.277491e-01,1.288020e-01,&
    1.298584e-01,1.309163e-01,1.319736e-01,1.330281e-01,1.340778e-01,&
    1.351207e-01,1.361546e-01,1.371775e-01,1.381873e-01,1.391820e-01,&
    1.401593e-01,1.411174e-01,1.420540e-01,1.429671e-01,1.438547e-01,&
    1.447146e-01,1.455449e-01,1.463433e-01,1.471078e-01,1.478364e-01,&
    1.485270e-01,1.491774e-01,1.497857e-01,1.503497e-01,1.508674e-01,&
    1.513367e-01,1.517554e-01,1.521216e-01,1.524332e-01,1.526880e-01,&
    1.528840e-01, &
    ! RRTMGP band 10 = RRTMG band 9 (12850-16000 cm-1)
    1.160842e-01,1.169118e-01,1.177697e-01,1.186556e-01,1.195676e-01,&
    1.205036e-01,1.214616e-01,1.224394e-01,1.234349e-01,1.244463e-01,&
    1.254712e-01,1.265078e-01,1.275539e-01,1.286075e-01,1.296664e-01,&
    1.307287e-01,1.317923e-01,1.328550e-01,1.339149e-01,1.349699e-01,&
    1.360179e-01,1.370567e-01,1.380845e-01,1.390991e-01,1.400984e-01,&
    1.410803e-01,1.420429e-01,1.429840e-01,1.439016e-01,1.447936e-01,&
    1.456579e-01,1.464925e-01,1.472953e-01,1.480642e-01,1.487972e-01,&
    1.494923e-01,1.501472e-01,1.507601e-01,1.513287e-01,1.518511e-01,&
    1.523252e-01,1.527489e-01,1.531201e-01,1.534368e-01,1.536969e-01,&
    1.538984e-01, &
    ! RRTMGP band 11 = RRTMG band 10 (16000-22650 cm-1)
    1.168725e-01,1.177088e-01,1.185747e-01,1.194680e-01,1.203867e-01,&
    1.213288e-01,1.222923e-01,1.232750e-01,1.242750e-01,1.252903e-01,&
    1.263187e-01,1.273583e-01,1.284069e-01,1.294626e-01,1.305233e-01,&
    1.315870e-01,1.326517e-01,1.337152e-01,1.347756e-01,1.358308e-01,&
    1.368788e-01,1.379175e-01,1.389449e-01,1.399590e-01,1.409577e-01,&
    1.419389e-01,1.429007e-01,1.438410e-01,1.447577e-01,1.456488e-01,&
    1.465123e-01,1.473461e-01,1.481483e-01,1.489166e-01,1.496492e-01,&
    1.503439e-01,1.509988e-01,1.516118e-01,1.521808e-01,1.527038e-01,&
    1.531788e-01,1.536037e-01,1.539764e-01,1.542951e-01,1.545575e-01,&
    1.547617e-01, &
    ! RRTMGP band 12 = RRTMG band 11 (22650-29000 cm-1)
    1.180509e-01,1.189025e-01,1.197820e-01,1.206875e-01,1.216171e-01,&
    1.225687e-01,1.235404e-01,1.245303e-01,1.255363e-01,1.265564e-01,&
    1.275888e-01,1.286313e-01,1.296821e-01,1.307392e-01,1.318006e-01,&
    1.328643e-01,1.339284e-01,1.349908e-01,1.360497e-01,1.371029e-01,&
    1.381486e-01,1.391848e-01,1.402095e-01,1.412208e-01,1.422165e-01,&
    1.431949e-01,1.441539e-01,1.450915e-01,1.460058e-01,1.468947e-01,&
    1.477564e-01,1.485888e-01,1.493900e-01,1.501580e-01,1.508907e-01,&
    1.515864e-01,1.522428e-01,1.528582e-01,1.534305e-01,1.539578e-01,&
    1.544380e-01,1.548692e-01,1.552494e-01,1.555767e-01,1.558490e-01,&
    1.560645e-01, &
    ! RRTMGP band 13 = RRTMG band 12 (29000-38000 cm-1)
    1.200480e-01,1.209267e-01,1.218304e-01,1.227575e-01,1.237059e-01,&
    1.246739e-01,1.256595e-01,1.266610e-01,1.276765e-01,1.287041e-01,&
    1.297420e-01,1.307883e-01,1.318412e-01,1.328988e-01,1.339593e-01,&
    1.350207e-01,1.360813e-01,1.371393e-01,1.381926e-01,1.392396e-01,&
    1.402783e-01,1.413069e-01,1.423235e-01,1.433263e-01,1.443134e-01,&
    1.452830e-01,1.462332e-01,1.471622e-01,1.480681e-01,1.489490e-01,&
    1.498032e-01,1.506286e-01,1.514236e-01,1.521863e-01,1.529147e-01,&
    1.536070e-01,1.542614e-01,1.548761e-01,1.554491e-01,1.559787e-01,&
    1.564629e-01,1.568999e-01,1.572879e-01,1.576249e-01,1.579093e-01,&
    1.581390e-01, &
    ! RRTMGP band 14 = RRTMG band 13 (38000-50000 cm-1)
    1.247813e-01,1.256496e-01,1.265417e-01,1.274560e-01,1.283905e-01,&
    1.293436e-01,1.303135e-01,1.312983e-01,1.322964e-01,1.333060e-01,&
    1.343252e-01,1.353523e-01,1.363855e-01,1.374231e-01,1.384632e-01,&
    1.395042e-01,1.405441e-01,1.415813e-01,1.426140e-01,1.436404e-01,&
    1.446587e-01,1.456672e-01,1.466640e-01,1.476475e-01,1.486157e-01,&
    1.495671e-01,1.504997e-01,1.514117e-01,1.523016e-01,1.531673e-01,&
    1.540073e-01,1.548197e-01,1.556026e-01,1.563545e-01,1.570734e-01,&
    1.577576e-01,1.584054e-01,1.590149e-01,1.595843e-01,1.601120e-01,&
    1.605962e-01,1.610349e-01,1.614266e-01,1.617693e-01,1.620614e-01,&
    1.623011e-01 &
    ],shape(fdlice3_rrtmgp)),kind=wp)
  ! ----------------------------------------------------------------

contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

  subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION: This version uses the MAPL\_GenericSetServices. This function sets
!   the Initialize and Finalize services, as well as allocating our instance of a
!   MAPL\_MetaComp and putting it in the gridded component (GC). Here we only need
!   to register the Run method with ESMF and register the state variable
!   specifications with MAPL. \newline
!

!EOP

!=============================================================================

! ErrLog Variables

    character(len=ESMF_MAXSTR) :: IAm
    character(len=ESMF_MAXSTR) :: COMP_NAME
    integer                    :: STATUS

! Local derived type aliases

    type (MAPL_MetaComp), pointer :: MAPL

! Locals

    integer :: RUN_DT
    integer :: MY_STEP
    integer :: ACCUMINT
    real    :: DT

    logical :: USE_RRTMGP, USE_RRTMG, USE_CHOU
    integer :: NUM_BANDS_SOLAR
    logical :: SOLAR_TO_OBIO

    type (ty_RRTMGP_state), pointer :: rrtmgp_state
    type (ty_RRTMGP_wrap)           :: wrap

!=============================================================================

    ! Get my name and set-up traceback handle
    call ESMF_GridCompGet(GC, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) // 'SetServices'

    ! save pointer to the wrapped RRTMGP internal state in the GC
    allocate(rrtmgp_state, __STAT__)
    wrap%ptr => rrtmgp_state
    call ESMF_UserCompSetInternalState(GC, 'RRTMGP_state', wrap, status)
    VERIFY_(status)

    ! Get my internal MAPL_Generic state
    call MAPL_GetObjectFromGC (GC, MAPL, __RC__)

    ! Get the intervals; "heartbeat" must exist
    call MAPL_GetResource (MAPL, DT, Label="RUN_DT:", __RC__)
    RUN_DT = nint(DT)

    ! Refresh interval defaults to heartbeat.
    call MAPL_GetResource (MAPL, DT, Label=trim(COMP_NAME)//"_DT:", default=DT, __RC__)
    MY_STEP = nint(DT)

    ! Averaging interval defaults to refresh interval.
    call MAPL_GetResource (MAPL, DT, Label=trim(COMP_NAME)//"Avrg:", default=DT, __RC__)
    ACCUMINT = nint(DT)

    ! Decide which radiation to use:
    ! Needed in SetServices because we Export a per-band flux and the
    !   number of bands differs between codes.
    !----------------------------------------------------------------------
    call choose_solar_scheme (MAPL, USE_RRTMGP, USE_RRTMG, USE_CHOU, __RC__)

    ! Set number of solar bands
    if (USE_RRTMGP) then
      NUM_BANDS_SOLAR = NB_RRTMGP
    else if (USE_RRTMG) then
      NUM_BANDS_SOLAR = NB_RRTMG
    else
      NUM_BANDS_SOLAR = NB_CHOU
    end if

    ! Decide if should make OBIO exports
    call MAPL_GetResource ( MAPL, DO_OBIO, Label="USE_OCEANOBIOGEOCHEM:",DEFAULT=0, RC=STATUS)
    VERIFY_(STATUS)
    
    SOLAR_TO_OBIO = (DO_OBIO/=0)

! Set the state variable specs.
! -----------------------------

!BOS

! !IMPORT STATE:

    call MAPL_AddImportSpec(GC,                                              &
       LONG_NAME  = 'air_pressure',                                          &
       UNITS      = 'Pa',                                                    &
       SHORT_NAME = 'PLE',                                                   &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                      &
       AVERAGING_INTERVAL = ACCUMINT,                                        &
       REFRESH_INTERVAL   = MY_STEP,                                   __RC__)

    call MAPL_AddImportSpec(GC,                                              &
       LONG_NAME          = 'surface_skin_temperature',                      &
       UNITS              = 'K',                                             &
       SHORT_NAME         = 'TS',                                            &
       DIMS               = MAPL_DimsHorzOnly,                               &
       VLOCATION          = MAPL_VLocationNone,                              &
       AVERAGING_INTERVAL = ACCUMINT,                                        &
       REFRESH_INTERVAL   = MY_STEP,                                   __RC__)

    call MAPL_AddImportSpec(GC,                                              &
       LONG_NAME          = 'methane_concentration',                         &
       UNITS              = 'pppv',                                          &
       SHORT_NAME         = 'CH4',                                           &
       DIMS               = MAPL_DimsHorzVert,                               &
       VLOCATION          = MAPL_VLocationCenter,                            &
       AVERAGING_INTERVAL = ACCUMINT,                                        &
       REFRESH_INTERVAL   = MY_STEP,                                   __RC__)

    call MAPL_AddImportSpec(GC,                                              &
       LONG_NAME          = 'nitrous_oxide_concentration',                   &
       UNITS              = 'pppv',                                          &
       SHORT_NAME         = 'N2O',                                           &
       DIMS               = MAPL_DimsHorzVert,                               &
       VLOCATION          = MAPL_VLocationCenter,                            &
       AVERAGING_INTERVAL = ACCUMINT,                                        &
       REFRESH_INTERVAL   = MY_STEP,                                   __RC__)

    call MAPL_AddImportSpec(GC,                                              &
       LONG_NAME  = 'air_temperature',                                       &
       UNITS      = 'K',                                                     &
       SHORT_NAME = 'T',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       AVERAGING_INTERVAL = ACCUMINT,                                        &
       REFRESH_INTERVAL   = MY_STEP,                                   __RC__)

    call MAPL_AddImportSpec(GC,                                              &
       LONG_NAME  = 'specific_humidity',                                     &
       UNITS      = 'kg kg-1',                                               &
       SHORT_NAME = 'QV',                                                    &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       AVERAGING_INTERVAL = ACCUMINT,                                        &
       REFRESH_INTERVAL   = MY_STEP,                                   __RC__)

    call MAPL_AddImportSpec(GC,                                              &
       LONG_NAME  = 'mass_fraction_of_cloud_liquid_water_in_air',            &
       UNITS      = 'kg kg-1',                                               &
       SHORT_NAME = 'QL',                                                    &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       AVERAGING_INTERVAL = ACCUMINT,                                        &
       REFRESH_INTERVAL   = MY_STEP,                                   __RC__)

    call MAPL_AddImportSpec(GC,                                              &
       LONG_NAME  = 'mass_fraction_of_cloud_ice_in_air',                     &
       UNITS      = 'kg kg-1',                                               &
       SHORT_NAME = 'QI',                                                    &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       AVERAGING_INTERVAL = ACCUMINT,                                        &
       REFRESH_INTERVAL   = MY_STEP,                                   __RC__)

    call MAPL_AddImportSpec(GC,                                              &
       LONG_NAME  = 'mass_fraction_of_rain_water_in_air',                    &
       UNITS      = 'kg kg-1',                                               &
       SHORT_NAME = 'QR',                                                    &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       AVERAGING_INTERVAL = ACCUMINT,                                        &
       REFRESH_INTERVAL   = MY_STEP,                                   __RC__)

    call MAPL_AddImportSpec(GC,                                              &
       LONG_NAME  = 'mass_fraction_of_snow_in_air',                          &
       UNITS      = 'kg kg-1',                                               &
       SHORT_NAME = 'QS',                                                    &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       AVERAGING_INTERVAL = ACCUMINT,                                        &
       REFRESH_INTERVAL   = MY_STEP,                                   __RC__)

    call MAPL_AddImportSpec(GC,                                              &
       LONG_NAME  = 'effective_radius_of_cloud_liquid_water_particles',      &
       UNITS      = 'm',                                                     &
       SHORT_NAME = 'RL',                                                    &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       AVERAGING_INTERVAL = ACCUMINT,                                        &
       REFRESH_INTERVAL   = MY_STEP,                                   __RC__)

    call MAPL_AddImportSpec(GC,                                              &
       LONG_NAME  = 'effective_radius_of_cloud_ice_particles',               &
       UNITS      = 'm',                                                     &
       SHORT_NAME = 'RI',                                                    &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       AVERAGING_INTERVAL = ACCUMINT,                                        &
       REFRESH_INTERVAL   = MY_STEP,                                   __RC__)

    call MAPL_AddImportSpec(GC,                                              &
       LONG_NAME  = 'effective_radius_of_rain_particles',                    &
       UNITS      = 'm',                                                     &
       SHORT_NAME = 'RR',                                                    &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       AVERAGING_INTERVAL = ACCUMINT,                                        &
       REFRESH_INTERVAL   = MY_STEP,                                   __RC__)

    call MAPL_AddImportSpec(GC,                                              &
       LONG_NAME  = 'effective_radius_of_snow_particles',                    &
       UNITS      = 'm',                                                     &
       SHORT_NAME = 'RS',                                                    &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       AVERAGING_INTERVAL = ACCUMINT,                                        &
       REFRESH_INTERVAL   = MY_STEP,                                   __RC__)

    call MAPL_AddImportSpec(GC,                                              &
       LONG_NAME  = 'odd-oxygen_volume_mixing_ratio',                        &
       UNITS      = 'mol mol-1',                                             &
       SHORT_NAME = 'OX',                                                    &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       AVERAGING_INTERVAL = ACCUMINT,                                        &
       REFRESH_INTERVAL   = MY_STEP,                                   __RC__)

    call MAPL_AddImportSpec(GC,                                              &
       LONG_NAME  = 'cloud_area_fraction',                                   &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'FCLD',                                                  &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       AVERAGING_INTERVAL = ACCUMINT,                                        &
       REFRESH_INTERVAL   = MY_STEP,                                   __RC__)

    call MAPL_AddImportSpec(GC,                                              &
       LONG_NAME  = 'aerosols',                                              &
       UNITS      = 'kg kg-1',                                               &
       SHORT_NAME = 'AERO',                                                  &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
       DATATYPE   = MAPL_StateItem,                                          &
       RESTART    = MAPL_RestartSkip,                                  __RC__)

    call MAPL_AddImportSpec(GC,                                              &
       LONG_NAME  = 'surface_albedo_for_visible_beam',                       &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'ALBVR',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddImportSpec(GC,                                              &
       LONG_NAME  = 'surface_albedo_for_visible_diffuse',                    &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'ALBVF',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddImportSpec(GC,                                              &
       LONG_NAME  = 'surface_albedo_for_near_infrared_beam',                 &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'ALBNR',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddImportSpec(GC,                                              &
       LONG_NAME  = 'surface_albedo_for_near_infrared_diffuse',              &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'ALBNF',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddImportSpec(GC,                                              &
       SHORT_NAME = 'PREF',                                                  &
       LONG_NAME  = 'reference_air_pressure',                                &
       UNITS      = 'Pa',                                                    &
       DIMS       = MAPL_DimsVertOnly,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                __RC__)


!  Solar does not have a "real" state. We keep an internal variable
!  for each variable produced by solar during the compute steps.
!  Versions of these, weighted by the appropriate TOA insolation,
!  are returned at each time step.

!  NB: Those INTERNALs with "FRIENDLYTO = trim(COMP_NAME)" effect a
!  MAPL feature where INTERNALs can be EXPORTed without an explicit
!  EXPORT of the same name and without any explicit INTERNAL to
!  EXPORT copying (such as we used to do with DRUVRN, etc., after
!  the second SORADCORE call in REFRESH). With this feature there is
!  only INTERNAL storage reserved, and the EXPORT, if requested, just
!  gets its value from that space.

!  !INTERNAL STATE:

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  ='normalized_net_downward_shortwave_flux_in_air',          &
       UNITS      ='1',                                                      &
       SHORT_NAME ='FSWN',                                                   &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                __RC__)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  ='normalized_net_downward_shortwave_flux_in_air_assuming_clear_sky',&
       UNITS      ='1',                                                      &
       SHORT_NAME ='FSCN',                                                   &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                __RC__)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  ='normalized_upward_shortwave_flux_in_air',                &
       UNITS      ='1',                                                      &
       SHORT_NAME ='FSWUN',                                                  &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                __RC__)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  ='normalized_upward_shortwave_flux_in_air_assuming_clear_sky',&
       UNITS      ='1',                                                      &
       SHORT_NAME ='FSCUN',                                                  &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                __RC__)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME      = 'normalized_net_surface_downward_shortwave_flux_per_band_in_air',&
       UNITS          = '1',                                                 &
       SHORT_NAME     = 'FSWBANDN',                                          &
       DIMS           = MAPL_DimsHorzOnly,                                   &
       UNGRIDDED_DIMS = (/ NUM_BANDS_SOLAR /),                               &
       VLOCATION      = MAPL_VLocationNone,                            __RC__)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'normalized_surface_downwelling_ultraviolet_beam_flux',  &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'DRUVRN',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
       FRIENDLYTO = trim(COMP_NAME),                                   __RC__)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'normalized_surface_downwelling_ultraviolet_diffuse_flux',&
       UNITS      = '1',                                                     &
       SHORT_NAME = 'DFUVRN',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
       FRIENDLYTO = trim(COMP_NAME),                                   __RC__)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'normalized_surface_downwelling_par_beam_flux',          &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'DRPARN',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
       FRIENDLYTO = trim(COMP_NAME),                                   __RC__)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'normalized_surface_downwelling_par_diffuse_flux',       &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'DFPARN',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
       FRIENDLYTO = trim(COMP_NAME),                                   __RC__)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'normalized_surface_downwelling_nearinfrared_beam_flux', &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'DRNIRN',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
       FRIENDLYTO = trim(COMP_NAME),                                   __RC__)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  = 'normalized_surface_downwelling_nearinfrared_diffuse_flux',&
       UNITS      = '1',                                                     &
       SHORT_NAME = 'DFNIRN',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
       FRIENDLYTO = trim(COMP_NAME),                                   __RC__)

    if (SOLAR_TO_OBIO) then

       call MAPL_AddInternalSpec(GC,                                         &
          LONG_NAME      = 'normalized_surface_downwelling_shortwave_beam_flux_per_band',&
          UNITS          = '1',                                              &
          SHORT_NAME     = 'DRBANDN',                                        &
          DIMS           = MAPL_DimsHorzOnly,                                &
          UNGRIDDED_DIMS = (/ NUM_BANDS_SOLAR /),                            &
          VLOCATION      = MAPL_VLocationNone,                         __RC__)

       call MAPL_AddInternalSpec(GC,                                         &
          LONG_NAME      = 'normalized_surface_downwelling_shortwave_diffuse_flux_per_band',&
          UNITS          = '1',                                              &
          SHORT_NAME     = 'DFBANDN',                                        &
          DIMS           = MAPL_DimsHorzOnly,                                &
          UNGRIDDED_DIMS = (/ NUM_BANDS_SOLAR /),                            &
          VLOCATION      = MAPL_VLocationNone,                         __RC__)

    end if

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  ='normalized_net_downward_shortwave_flux_in_air_assuming_no_aerosol',&
       UNITS      ='1',                                                      &
       SHORT_NAME ='FSWNAN',                                                 &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                __RC__)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  ='normalized_net_downward_shortwave_flux_in_air_assuming_clear_sky_and_no_aerosol',&
       UNITS      ='1',                                                      &
       SHORT_NAME ='FSCNAN',                                                 &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                __RC__)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  ='normalized_upward_shortwave_flux_in_air_assuming_no_aerosol',&
       UNITS      ='1',                                                      &
       SHORT_NAME ='FSWUNAN',                                                &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                __RC__)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME  ='normalized_upward_shortwave_flux_in_air_assuming_clear_sky_and_no_aerosol',&
       UNITS      ='1',                                                      &
       SHORT_NAME ='FSCUNAN',                                                &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                __RC__)

    call MAPL_AddInternalSpec(GC,                                            &
       LONG_NAME      = 'normalized_net_surface_downward_shortwave_flux_per_band_in_air_assuming_no_aerosol',&
       UNITS          = '1',                                                 &
       SHORT_NAME     = 'FSWBANDNAN',                                        &
       DIMS           = MAPL_DimsHorzOnly,                                   &
       UNGRIDDED_DIMS = (/ NUM_BANDS_SOLAR /),                               &
       VLOCATION      = MAPL_VLocationNone,                            __RC__)

!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  NB: The following INTERNALs are really EXPORTs. As of 5/2022 MAPL only re-
!  loads IMPORTs and INTERNALs at the beginning of the replay Corrector phase,
!  not EXPORTs, with the consequence that intermittant EXPORTs (such as these
!  solar REFRESH exports) which were assumed to persist will actually have
!  values from the end of the Predictor phase, which is 3 hours ahead.
!    The fix below uses the combined INTERNAL/EXPORT feature of MAPL mentioned
!  above (with the required FRIENDLYTO). Since the EXPORT becomes an INTERNAL,
!  it IS reloaded and time-synced correctly. It is not a perfect solution,
!  since now EXPORTS that are NOT requested still take up INTERNAL storage.
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! This is the flux-weighted and time-averaged value over the REFRESH interval.
    ! It should be used for sub-sampling or analysing REFRESH diagnostics, such
    ! as CLDxxSW.
    call MAPL_AddInternalSpec(GC,                                            &
       SHORT_NAME = 'COSZSW',                                                &
       LONG_NAME  = 'cosine_of_the_solar_zenith_angle_of_Solar_REFRESH',     &
       UNITS      = '1',                                                     &
       DEFAULT    = MAPL_UNDEF,                                              &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
       FRIENDLYTO = trim(COMP_NAME),                                   __RC__)

    ! Note: the four CLDxxSW diagnostics below represent super-layer cloud
    ! fractions based on the subcolumn cloud generation called in RRTMG[P] SW.
    ! They are sunlit only fields and generated only at the SW REFRESH frequency,
    ! NOT at the heartbeat. As such, they are useful for diagnostic comparisons
    ! with the CLDxx set above. But they should NOT be used to subsample fields
    ! that are produced on the model heartbeat (e.g. subsampling for cloud
    ! presence). Note, also, that when comparing CLDxxSW with CLDxx, it is better
    ! to subsample both with COSZSW >= cmin, (e.g., 0.25). This COSZSW is a
    ! REFRESH-frequency version of MCOSZ and, as such, is most appropriate for
    ! subsampling REFRESH-frequency fields like CLDxxSW. Of course, you can also
    ! subsample CLDxx with COSZSW since CLDxx are global. By sampling both
    ! CLDxxSW and CLDxx with COSZSW you get a fair apples-to-apples comparison
    ! between the two.

    call MAPL_AddInternalSpec(GC,                                            &
       SHORT_NAME = 'CLDTTSW',                                               &
       LONG_NAME  = 'total_cloud_area_fraction_RRTMG_P_SW_REFRESH',          &
       UNITS      = '1',                                                     &
       DEFAULT    = MAPL_UNDEF,                                              &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
       FRIENDLYTO = trim(COMP_NAME),                                   __RC__)

    call MAPL_AddInternalSpec(GC,                                            &
       SHORT_NAME = 'CLDHISW',                                               &
       LONG_NAME  = 'high-level_cloud_area_fraction_RRTMG_P_SW_REFRESH',     &
       UNITS      = '1',                                                     &
       DEFAULT    = MAPL_UNDEF,                                              &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
       FRIENDLYTO = trim(COMP_NAME),                                   __RC__)

    call MAPL_AddInternalSpec(GC,                                            &
       SHORT_NAME = 'CLDMDSW',                                               &
       LONG_NAME  = 'mid-level_cloud_area_fraction_RRTMG_P_SW_REFRESH',      &
       UNITS      = '1',                                                     &
       DEFAULT    = MAPL_UNDEF,                                              &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
       FRIENDLYTO = trim(COMP_NAME),                                   __RC__)

    call MAPL_AddInternalSpec(GC,                                            &
       SHORT_NAME = 'CLDLOSW',                                               &
       LONG_NAME  = 'low-level_cloud_area_fraction_RRTMG_P_SW_REFRESH',      &
       UNITS      = '1',                                                     &
       DEFAULT    = MAPL_UNDEF,                                              &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
       FRIENDLYTO = trim(COMP_NAME),                                   __RC__)

    ! Note: The following TAUxxPAR and COTxxPAR are REFRESH-frequency fields.
    ! As such, all the important provisos given in the comment on CLDxxSW above
    ! apply to these fields as well. Please read those provisos. Their advantage
    ! is that they use the subcolumn cloud generation called in RRTMG[P].

#ifdef SOLAR_RADVAL
    ! TAUxxPAR are ZERO for clear super-layers, an anti-pattern for *in-cloud* optical
    ! thicknesses, and deprecated. They are currently included under the SOLAR_RADVAL flag,
    ! which is generally reserved for developer usage. They may later be removed completely.

    call MAPL_AddInternalSpec(GC,                                                      &
       SHORT_NAME = 'TAULOPAR',                                                        &
       LONG_NAME  = 'in_cloud_optical_thickness_of_low_clouds_RRTMG_P_PAR_REFRESH__deprecated', &
       UNITS      = '1' ,                                                              &
       DEFAULT    = MAPL_UNDEF,                                                        &
       DIMS       = MAPL_DimsHorzOnly,                                                 &
       VLOCATION  = MAPL_VLocationNone,                                                &
       FRIENDLYTO = trim(COMP_NAME),                                             __RC__)

    call MAPL_AddInternalSpec(GC,                                                      &
       SHORT_NAME = 'TAUMDPAR',                                                        &
       LONG_NAME  = 'in_cloud_optical_thickness_of_middle_clouds_RRTMG_P_PAR_REFRESH__deprecated', &
       UNITS      = '1' ,                                                              &
       DEFAULT    = MAPL_UNDEF,                                                        &
       DIMS       = MAPL_DimsHorzOnly,                                                 &
       VLOCATION  = MAPL_VLocationNone,                                                &
       FRIENDLYTO = trim(COMP_NAME),                                             __RC__)

    call MAPL_AddInternalSpec(GC,                                                      &
       SHORT_NAME = 'TAUHIPAR',                                                        &
       LONG_NAME  = 'in_cloud_optical_thickness_of_high_clouds_RRTMG_P_PAR_REFRESH__deprecated', &
       UNITS      = '1' ,                                                              &
       DEFAULT    = MAPL_UNDEF,                                                        &
       DIMS       = MAPL_DimsHorzOnly,                                                 &
       VLOCATION  = MAPL_VLocationNone,                                                &
       FRIENDLYTO = trim(COMP_NAME),                                             __RC__)

    call MAPL_AddInternalSpec(GC,                                                      &
       SHORT_NAME = 'TAUTTPAR',                                                        &
       LONG_NAME  = 'in_cloud_optical_thickness_of_all_clouds_RRTMG_P_PAR_REFRESH__deprecated', &
       UNITS      = '1' ,                                                              &
       DEFAULT    = MAPL_UNDEF,                                                        &
       DIMS       = MAPL_DimsHorzOnly,                                                 &
       VLOCATION  = MAPL_VLocationNone,                                                &
       FRIENDLYTO = trim(COMP_NAME),                                             __RC__)
#endif

    ! These COTxxPAR are UNDEF for clear super-layers, whereas TAUxxPAR are ZERO.
    ! As such, the COTxxPAR are a better in-cloud diagnostic.

    call MAPL_AddInternalSpec(GC,                                                      &
       SHORT_NAME = 'COTLOPAR',                                                        &
       LONG_NAME  = 'in_cloud_optical_thickness_of_low_clouds_RRTMG_P_PAR_REFRESH_clrundef', &
       UNITS      = '1' ,                                                              &
       DEFAULT    = MAPL_UNDEF,                                                        &
       DIMS       = MAPL_DimsHorzOnly,                                                 &
       VLOCATION  = MAPL_VLocationNone,                                                &
       FRIENDLYTO = trim(COMP_NAME),                                             __RC__)

    call MAPL_AddInternalSpec(GC,                                                      &
       SHORT_NAME = 'COTMDPAR',                                                        &
       LONG_NAME  = 'in_cloud_optical_thickness_of_middle_clouds_RRTMG_P_PAR_REFRESH_clrundef', &
       UNITS      = '1' ,                                                              &
       DEFAULT    = MAPL_UNDEF,                                                        &
       DIMS       = MAPL_DimsHorzOnly,                                                 &
       VLOCATION  = MAPL_VLocationNone,                                                &
       FRIENDLYTO = trim(COMP_NAME),                                             __RC__)

    call MAPL_AddInternalSpec(GC,                                                      &
       SHORT_NAME = 'COTHIPAR',                                                        &
       LONG_NAME  = 'in_cloud_optical_thickness_of_high_clouds_RRTMG_P_PAR_REFRESH_clrundef', &
       UNITS      = '1' ,                                                              &
       DEFAULT    = MAPL_UNDEF,                                                        &
       DIMS       = MAPL_DimsHorzOnly,                                                 &
       VLOCATION  = MAPL_VLocationNone,                                                &
       FRIENDLYTO = trim(COMP_NAME),                                             __RC__)

    call MAPL_AddInternalSpec(GC,                                                      &
       SHORT_NAME = 'COTTTPAR',                                                        &
       LONG_NAME  = 'in_cloud_optical_thickness_of_all_clouds_RRTMG_P_PAR_REFRESH_clrundef', &
       UNITS      = '1' ,                                                              &
       DEFAULT    = MAPL_UNDEF,                                                        &
       DIMS       = MAPL_DimsHorzOnly,                                                 &
       VLOCATION  = MAPL_VLocationNone,                                                &
       FRIENDLYTO = trim(COMP_NAME),                                             __RC__)

    ! For COT[DEN|NUM]xxPAR see comments under COT[DEN|NUM]xx.

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTDENLOPAR',                                                                 &
       LONG_NAME  = 'in_cloud_optical_thickness_of_low_clouds_RRTMG_P_PAR_REFRESH_denominator',    &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTDENMDPAR',                                                                 &
       LONG_NAME  = 'in_cloud_optical_thickness_of_middle_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTDENHIPAR',                                                                 &
       LONG_NAME  = 'in_cloud_optical_thickness_of_high_clouds_RRTMG_P_PAR_REFRESH_denominator',   &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTDENTTPAR',                                                                 &
       LONG_NAME  = 'in_cloud_optical_thickness_of_all_clouds_RRTMG_P_PAR_REFRESH_denominator',    &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTNUMLOPAR',                                                                 &
       LONG_NAME  = 'in_cloud_optical_thickness_of_low_clouds_RRTMG_P_PAR_REFRESH_numerator',      &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTNUMMDPAR',                                                                 &
       LONG_NAME  = 'in_cloud_optical_thickness_of_middle_clouds_RRTMG_P_PAR_REFRESH_numerator',   &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTNUMHIPAR',                                                                 &
       LONG_NAME  = 'in_cloud_optical_thickness_of_high_clouds_RRTMG_P_PAR_REFRESH_numerator',     &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTNUMTTPAR',                                                                 &
       LONG_NAME  = 'in_cloud_optical_thickness_of_all_clouds_RRTMG_P_PAR_REFRESH_numerator',      &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

#ifdef SOLAR_RADVAL

    ! COTDS[DEN|NUM]xxPAR are like COT[DEN|NUM]xxPAR but delta-scaled.

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTDSDENLOPAR',                                                               &
       LONG_NAME  = 'in_cloud_optical_thickness_delta_scaled_of_low_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTDSDENMDPAR',                                                               &
       LONG_NAME  = 'in_cloud_optical_thickness_delta_scaled_of_middle_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTDSDENHIPAR',                                                               &
       LONG_NAME  = 'in_cloud_optical_thickness_delta_scaled_of_high_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTDSDENTTPAR',                                                               &
       LONG_NAME  = 'in_cloud_optical_thickness_delta_scaled_of_all_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTDSNUMLOPAR',                                                               &
       LONG_NAME  = 'in_cloud_optical_thickness_delta_scaled_of_low_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTDSNUMMDPAR',                                                               &
       LONG_NAME  = 'in_cloud_optical_thickness_delta_scaled_of_middle_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTDSNUMHIPAR',                                                               &
       LONG_NAME  = 'in_cloud_optical_thickness_delta_scaled_of_high_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTDSNUMTTPAR',                                                               &
       LONG_NAME  = 'in_cloud_optical_thickness_delta_scaled_of_all_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    ! ditto for liquid clouds only

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTLDENLOPAR',                                                                &
       LONG_NAME  = 'in_cloud_liquid_optical_thickness_of_low_clouds_RRTMG_P_PAR_REFRESH_denominator',    &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTLDENMDPAR',                                                                &
       LONG_NAME  = 'in_cloud_liquid_optical_thickness_of_middle_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTLDENHIPAR',                                                                &
       LONG_NAME  = 'in_cloud_liquid_optical_thickness_of_high_clouds_RRTMG_P_PAR_REFRESH_denominator',   &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTLDENTTPAR',                                                                &
       LONG_NAME  = 'in_cloud_liquid_optical_thickness_of_all_clouds_RRTMG_P_PAR_REFRESH_denominator',    &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTLNUMLOPAR',                                                                &
       LONG_NAME  = 'in_cloud_liquid_optical_thickness_of_low_clouds_RRTMG_P_PAR_REFRESH_numerator',      &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTLNUMMDPAR',                                                                &
       LONG_NAME  = 'in_cloud_liquid_optical_thickness_of_middle_clouds_RRTMG_P_PAR_REFRESH_numerator',   &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTLNUMHIPAR',                                                                &
       LONG_NAME  = 'in_cloud_liquid_optical_thickness_of_high_clouds_RRTMG_P_PAR_REFRESH_numerator',     &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTLNUMTTPAR',                                                                &
       LONG_NAME  = 'in_cloud_liquid_optical_thickness_of_all_clouds_RRTMG_P_PAR_REFRESH_numerator',      &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTLDSDENLOPAR',                                                              &
       LONG_NAME  = 'in_cloud_liquid_optical_thickness_delta_scaled_of_low_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTLDSDENMDPAR',                                                              &
       LONG_NAME  = 'in_cloud_liquid_optical_thickness_delta_scaled_of_middle_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTLDSDENHIPAR',                                                              &
       LONG_NAME  = 'in_cloud_liquid_optical_thickness_delta_scaled_of_high_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTLDSDENTTPAR',                                                              &
       LONG_NAME  = 'in_cloud_liquid_optical_thickness_delta_scaled_of_all_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTLDSNUMLOPAR',                                                              &
       LONG_NAME  = 'in_cloud_liquid_optical_thickness_delta_scaled_of_low_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTLDSNUMMDPAR',                                                              &
       LONG_NAME  = 'in_cloud_liquid_optical_thickness_delta_scaled_of_middle_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTLDSNUMHIPAR',                                                              &
       LONG_NAME  = 'in_cloud_liquid_optical_thickness_delta_scaled_of_high_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTLDSNUMTTPAR',                                                              &
       LONG_NAME  = 'in_cloud_liquid_optical_thickness_delta_scaled_of_all_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    ! ditto for ice clouds only

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTIDENLOPAR',                                                                &
       LONG_NAME  = 'in_cloud_ice_optical_thickness_of_low_clouds_RRTMG_P_PAR_REFRESH_denominator',    &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTIDENMDPAR',                                                                &
       LONG_NAME  = 'in_cloud_ice_optical_thickness_of_middle_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTIDENHIPAR',                                                                &
       LONG_NAME  = 'in_cloud_ice_optical_thickness_of_high_clouds_RRTMG_P_PAR_REFRESH_denominator',   &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTIDENTTPAR',                                                                &
       LONG_NAME  = 'in_cloud_ice_optical_thickness_of_all_clouds_RRTMG_P_PAR_REFRESH_denominator',    &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTINUMLOPAR',                                                                &
       LONG_NAME  = 'in_cloud_ice_optical_thickness_of_low_clouds_RRTMG_P_PAR_REFRESH_numerator',      &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTINUMMDPAR',                                                                &
       LONG_NAME  = 'in_cloud_ice_optical_thickness_of_middle_clouds_RRTMG_P_PAR_REFRESH_numerator',   &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTINUMHIPAR',                                                                &
       LONG_NAME  = 'in_cloud_ice_optical_thickness_of_high_clouds_RRTMG_P_PAR_REFRESH_numerator',     &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTINUMTTPAR',                                                                &
       LONG_NAME  = 'in_cloud_ice_optical_thickness_of_all_clouds_RRTMG_P_PAR_REFRESH_numerator',      &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTIDSDENLOPAR',                                                              &
       LONG_NAME  = 'in_cloud_ice_optical_thickness_delta_scaled_of_low_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTIDSDENMDPAR',                                                              &
       LONG_NAME  = 'in_cloud_ice_optical_thickness_delta_scaled_of_middle_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTIDSDENHIPAR',                                                              &
       LONG_NAME  = 'in_cloud_ice_optical_thickness_delta_scaled_of_high_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTIDSDENTTPAR',                                                              &
       LONG_NAME  = 'in_cloud_ice_optical_thickness_delta_scaled_of_all_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTIDSNUMLOPAR',                                                              &
       LONG_NAME  = 'in_cloud_ice_optical_thickness_delta_scaled_of_low_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTIDSNUMMDPAR',                                                              &
       LONG_NAME  = 'in_cloud_ice_optical_thickness_delta_scaled_of_middle_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTIDSNUMHIPAR',                                                              &
       LONG_NAME  = 'in_cloud_ice_optical_thickness_delta_scaled_of_high_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'COTIDSNUMTTPAR',                                                              &
       LONG_NAME  = 'in_cloud_ice_optical_thickness_delta_scaled_of_all_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    ! super-layerized phase-split cloud SSA and ASM 

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'SSALDENLOPAR',                                                                &
       LONG_NAME  = 'in_cloud_liquid_single_scattering_albedo_of_low_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'SSALNUMLOPAR',                                                                &
       LONG_NAME  = 'in_cloud_liquid_single_scattering_albedo_of_low_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'SSAIDENLOPAR',                                                                &
       LONG_NAME  = 'in_cloud_ice_single_scattering_albedo_of_low_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'SSAINUMLOPAR',                                                                &
       LONG_NAME  = 'in_cloud_ice_single_scattering_albedo_of_low_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'ASMLDENLOPAR',                                                                &
       LONG_NAME  = 'in_cloud_liquid_asymmetry_parameter_of_low_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'ASMLNUMLOPAR',                                                                &
       LONG_NAME  = 'in_cloud_liquid_asymmetry_parameter_of_low_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'ASMIDENLOPAR',                                                                &
       LONG_NAME  = 'in_cloud_ice_asymmetry_parameter_of_low_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'ASMINUMLOPAR',                                                                &
       LONG_NAME  = 'in_cloud_ice_asymmetry_parameter_of_low_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'SSALDENMDPAR',                                                                &
       LONG_NAME  = 'in_cloud_liquid_single_scattering_albedo_of_middle_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'SSALNUMMDPAR',                                                                &
       LONG_NAME  = 'in_cloud_liquid_single_scattering_albedo_of_middle_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'SSAIDENMDPAR',                                                                &
       LONG_NAME  = 'in_cloud_ice_single_scattering_albedo_of_middle_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'SSAINUMMDPAR',                                                                &
       LONG_NAME  = 'in_cloud_ice_single_scattering_albedo_of_middle_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'ASMLDENMDPAR',                                                                &
       LONG_NAME  = 'in_cloud_liquid_asymmetry_parameter_of_middle_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'ASMLNUMMDPAR',                                                                &
       LONG_NAME  = 'in_cloud_liquid_asymmetry_parameter_of_middle_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'ASMIDENMDPAR',                                                                &
       LONG_NAME  = 'in_cloud_ice_asymmetry_parameter_of_middle_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'ASMINUMMDPAR',                                                                &
       LONG_NAME  = 'in_cloud_ice_asymmetry_parameter_of_middle_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'SSALDENHIPAR',                                                                &
       LONG_NAME  = 'in_cloud_liquid_single_scattering_albedo_of_high_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'SSALNUMHIPAR',                                                                &
       LONG_NAME  = 'in_cloud_liquid_single_scattering_albedo_of_high_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'SSAIDENHIPAR',                                                                &
       LONG_NAME  = 'in_cloud_ice_single_scattering_albedo_of_high_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'SSAINUMHIPAR',                                                                &
       LONG_NAME  = 'in_cloud_ice_single_scattering_albedo_of_high_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'ASMLDENHIPAR',                                                                &
       LONG_NAME  = 'in_cloud_liquid_asymmetry_parameter_of_high_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'ASMLNUMHIPAR',                                                                &
       LONG_NAME  = 'in_cloud_liquid_asymmetry_parameter_of_high_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'ASMIDENHIPAR',                                                                &
       LONG_NAME  = 'in_cloud_ice_asymmetry_parameter_of_high_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'ASMINUMHIPAR',                                                                &
       LONG_NAME  = 'in_cloud_ice_asymmetry_parameter_of_high_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'SSALDENTTPAR',                                                                &
       LONG_NAME  = 'in_cloud_liquid_single_scattering_albedo_of_all_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'SSALNUMTTPAR',                                                                &
       LONG_NAME  = 'in_cloud_liquid_single_scattering_albedo_of_all_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'SSAIDENTTPAR',                                                                &
       LONG_NAME  = 'in_cloud_ice_single_scattering_albedo_of_all_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'SSAINUMTTPAR',                                                                &
       LONG_NAME  = 'in_cloud_ice_single_scattering_albedo_of_all_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'ASMLDENTTPAR',                                                                &
       LONG_NAME  = 'in_cloud_liquid_asymmetry_parameter_of_all_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'ASMLNUMTTPAR',                                                                &
       LONG_NAME  = 'in_cloud_liquid_asymmetry_parameter_of_all_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'ASMIDENTTPAR',                                                                &
       LONG_NAME  = 'in_cloud_ice_asymmetry_parameter_of_all_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'ASMINUMTTPAR',                                                                &
       LONG_NAME  = 'in_cloud_ice_asymmetry_parameter_of_all_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'SSALDSDENLOPAR',                                                              &
       LONG_NAME  = 'in_cloud_liquid_single_scattering_albedo_delta_scaled_of_low_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'SSALDSNUMLOPAR',                                                              &
       LONG_NAME  = 'in_cloud_liquid_single_scattering_albedo_delta_scaled_of_low_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'SSAIDSDENLOPAR',                                                              &
       LONG_NAME  = 'in_cloud_ice_single_scattering_albedo_delta_scaled_of_low_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'SSAIDSNUMLOPAR',                                                              &
       LONG_NAME  = 'in_cloud_ice_single_scattering_albedo_delta_scaled_of_low_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'ASMLDSDENLOPAR',                                                              &
       LONG_NAME  = 'in_cloud_liquid_asymmetry_parameter_delta_scaled_of_low_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'ASMLDSNUMLOPAR',                                                              &
       LONG_NAME  = 'in_cloud_liquid_asymmetry_parameter_delta_scaled_of_low_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'ASMIDSDENLOPAR',                                                              &
       LONG_NAME  = 'in_cloud_ice_asymmetry_parameter_delta_scaled_of_low_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'ASMIDSNUMLOPAR',                                                              &
       LONG_NAME  = 'in_cloud_ice_asymmetry_parameter_delta_scaled_of_low_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'SSALDSDENMDPAR',                                                              &
       LONG_NAME  = 'in_cloud_liquid_single_scattering_albedo_delta_scaled_of_middle_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'SSALDSNUMMDPAR',                                                              &
       LONG_NAME  = 'in_cloud_liquid_single_scattering_albedo_delta_scaled_of_middle_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'SSAIDSDENMDPAR',                                                              &
       LONG_NAME  = 'in_cloud_ice_single_scattering_albedo_delta_scaled_of_middle_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'SSAIDSNUMMDPAR',                                                              &
       LONG_NAME  = 'in_cloud_ice_single_scattering_albedo_delta_scaled_of_middle_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'ASMLDSDENMDPAR',                                                              &
       LONG_NAME  = 'in_cloud_liquid_asymmetry_parameter_delta_scaled_of_middle_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'ASMLDSNUMMDPAR',                                                              &
       LONG_NAME  = 'in_cloud_liquid_asymmetry_parameter_delta_scaled_of_middle_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'ASMIDSDENMDPAR',                                                              &
       LONG_NAME  = 'in_cloud_ice_asymmetry_parameter_delta_scaled_of_middle_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'ASMIDSNUMMDPAR',                                                              &
       LONG_NAME  = 'in_cloud_ice_asymmetry_parameter_delta_scaled_of_middle_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'SSALDSDENHIPAR',                                                              &
       LONG_NAME  = 'in_cloud_liquid_single_scattering_albedo_delta_scaled_of_high_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'SSALDSNUMHIPAR',                                                              &
       LONG_NAME  = 'in_cloud_liquid_single_scattering_albedo_delta_scaled_of_high_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'SSAIDSDENHIPAR',                                                              &
       LONG_NAME  = 'in_cloud_ice_single_scattering_albedo_delta_scaled_of_high_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'SSAIDSNUMHIPAR',                                                              &
       LONG_NAME  = 'in_cloud_ice_single_scattering_albedo_delta_scaled_of_high_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'ASMLDSDENHIPAR',                                                              &
       LONG_NAME  = 'in_cloud_liquid_asymmetry_parameter_delta_scaled_of_high_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'ASMLDSNUMHIPAR',                                                              &
       LONG_NAME  = 'in_cloud_liquid_asymmetry_parameter_delta_scaled_of_high_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'ASMIDSDENHIPAR',                                                              &
       LONG_NAME  = 'in_cloud_ice_asymmetry_parameter_delta_scaled_of_high_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'ASMIDSNUMHIPAR',                                                              &
       LONG_NAME  = 'in_cloud_ice_asymmetry_parameter_delta_scaled_of_high_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'SSALDSDENTTPAR',                                                              &
       LONG_NAME  = 'in_cloud_liquid_single_scattering_albedo_delta_scaled_of_all_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'SSALDSNUMTTPAR',                                                              &
       LONG_NAME  = 'in_cloud_liquid_single_scattering_albedo_delta_scaled_of_all_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'SSAIDSDENTTPAR',                                                              &
       LONG_NAME  = 'in_cloud_ice_single_scattering_albedo_delta_scaled_of_all_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'SSAIDSNUMTTPAR',                                                              &
       LONG_NAME  = 'in_cloud_ice_single_scattering_albedo_delta_scaled_of_all_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'ASMLDSDENTTPAR',                                                              &
       LONG_NAME  = 'in_cloud_liquid_asymmetry_parameter_delta_scaled_of_all_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'ASMLDSNUMTTPAR',                                                              &
       LONG_NAME  = 'in_cloud_liquid_asymmetry_parameter_delta_scaled_of_all_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'ASMIDSDENTTPAR',                                                              &
       LONG_NAME  = 'in_cloud_ice_asymmetry_parameter_delta_scaled_of_all_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'ASMIDSNUMTTPAR',                                                              &
       LONG_NAME  = 'in_cloud_ice_asymmetry_parameter_delta_scaled_of_all_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    ! super-layerized phase-split cloud forward-scattering fraction

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'FORLDENLOPAR',                                                                &
       LONG_NAME  = 'in_cloud_liquid_forward_scattering_fraction_of_low_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'FORLNUMLOPAR',                                                                &
       LONG_NAME  = 'in_cloud_liquid_forward_scattering_fraction_of_low_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'FORIDENLOPAR',                                                                &
       LONG_NAME  = 'in_cloud_ice_forward_scattering_fraction_of_low_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'FORINUMLOPAR',                                                                &
       LONG_NAME  = 'in_cloud_ice_forward_scattering_fraction_of_low_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'FORLDENMDPAR',                                                                &
       LONG_NAME  = 'in_cloud_liquid_forward_scattering_fraction_of_mid_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'FORLNUMMDPAR',                                                                &
       LONG_NAME  = 'in_cloud_liquid_forward_scattering_fraction_of_mid_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'FORIDENMDPAR',                                                                &
       LONG_NAME  = 'in_cloud_ice_forward_scattering_fraction_of_mid_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'FORINUMMDPAR',                                                                &
       LONG_NAME  = 'in_cloud_ice_forward_scattering_fraction_of_mid_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'FORLDENHIPAR',                                                                &
       LONG_NAME  = 'in_cloud_liquid_forward_scattering_fraction_of_high_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'FORLNUMHIPAR',                                                                &
       LONG_NAME  = 'in_cloud_liquid_forward_scattering_fraction_of_high_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'FORIDENHIPAR',                                                                &
       LONG_NAME  = 'in_cloud_ice_forward_scattering_fraction_of_high_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'FORINUMHIPAR',                                                                &
       LONG_NAME  = 'in_cloud_ice_forward_scattering_fraction_of_high_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'FORLDENTTPAR',                                                                &
       LONG_NAME  = 'in_cloud_liquid_forward_scattering_fraction_of_all_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'FORLNUMTTPAR',                                                                &
       LONG_NAME  = 'in_cloud_liquid_forward_scattering_fraction_of_all_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'FORIDENTTPAR',                                                                &
       LONG_NAME  = 'in_cloud_ice_forward_scattering_fraction_of_all_clouds_RRTMG_P_PAR_REFRESH_denominator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)

    call MAPL_AddInternalSpec(GC,                                                                  &
       SHORT_NAME = 'FORINUMTTPAR',                                                                &
       LONG_NAME  = 'in_cloud_ice_forward_scattering_fraction_of_all_clouds_RRTMG_P_PAR_REFRESH_numerator', &
       UNITS      = '1' ,                                                                          &
       DIMS       = MAPL_DimsHorzOnly,                                                             &
       VLOCATION  = MAPL_VLocationNone,                                                            &
       FRIENDLYTO = trim(COMP_NAME),                                                         __RC__)
#endif

!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  END of EXPORTs masquerading as INTERNALs
!  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!  !EXPORT STATE:

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  ='net_downward_shortwave_flux_in_air',                     &
       UNITS      ='W m-2',                                                  &
       SHORT_NAME ='FSW',                                                    &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  ='net_downward_shortwave_flux_in_air_assuming_clear_sky',  &
       UNITS      ='W m-2',                                                  &
       SHORT_NAME ='FSC',                                                    &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  ='net_downward_shortwave_flux_in_air_assuming_no_aerosol', &
       UNITS      ='W m-2',                                                  &
       SHORT_NAME ='FSWNA',                                                  &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  ='net_downward_shortwave_flux_in_air_assuming_clear_sky_and_no_aerosol',&
       UNITS      ='W m-2',                                                  &
       SHORT_NAME ='FSCNA',                                                  &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  ='downward_shortwave_flux_in_air',                         &
       UNITS      ='W m-2',                                                  &
       SHORT_NAME ='FSWD',                                                   &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  ='downward_shortwave_flux_in_air_assuming_clear_sky',      &
       UNITS      ='W m-2',                                                  &
       SHORT_NAME ='FSCD',                                                   &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  ='downward_shortwave_flux_in_air_assuming_no_aerosol',     &
       UNITS      ='W m-2',                                                  &
       SHORT_NAME ='FSWDNA',                                                 &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  ='downward_shortwave_flux_in_air_assuming_clear_sky_and_no_aerosol',&
       UNITS      ='W m-2',                                                  &
       SHORT_NAME ='FSCDNA',                                                 &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  ='upward_shortwave_flux_in_air',                           &
       UNITS      ='W m-2',                                                  &
       SHORT_NAME ='FSWU',                                                   &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  ='upward_shortwave_flux_in_air_assuming_clear_sky',        &
       UNITS      ='W m-2',                                                  &
       SHORT_NAME ='FSCU',                                                   &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  ='upward_shortwave_flux_in_air_assuming_no_aerosol',       &
       UNITS      ='W m-2',                                                  &
       SHORT_NAME ='FSWUNA',                                                 &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  ='upward_shortwave_flux_in_air_assuming_clear_sky_and_no_aerosol',&
       UNITS      ='W m-2',                                                  &
       SHORT_NAME ='FSCUNA',                                                 &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationEdge,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME      = 'net_surface_downward_shortwave_flux_per_band_in_air',&
       UNITS          = 'W m-2',                                             &
       SHORT_NAME     = 'FSWBAND',                                           &
       DIMS           = MAPL_DimsHorzOnly,                                   &
       UNGRIDDED_DIMS = (/ NUM_BANDS_SOLAR /),                               &
       VLOCATION      = MAPL_VLocationNone,                            __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME      = 'net_surface_downward_shortwave_flux_per_band_in_air_assuming_no_aerosol',&
       UNITS          = 'W m-2',                                             &
       SHORT_NAME     = 'FSWBANDNA',                                         &
       DIMS           = MAPL_DimsHorzOnly,                                   &
       UNGRIDDED_DIMS = (/ NUM_BANDS_SOLAR /),                               &
       VLOCATION      = MAPL_VLocationNone,                            __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'surface_downwelling_ultraviolet_beam_normal_flux',      &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'DRNUVR',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'surface_downwelling_par_beam_normal_flux',              &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'DRNPAR',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'surface_downwelling_nearinfrared_beam_normal_flux',     &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'DRNNIR',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'surface_downwelling_ultraviolet_beam_flux',             &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'DRUVR',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'surface_downwelling_ultraviolet_diffuse_flux',          &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'DFUVR',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'surface_downwelling_par_beam_flux',                     &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'DRPAR',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                               __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'surface_downwelling_par_diffuse_flux',                  &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'DFPAR',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'surface_downwelling_nearinfrared_beam_flux',            &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'DRNIR',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'surface_downwelling_nearinfrared_diffuse_flux',         &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'DFNIR',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    if (SOLAR_TO_OBIO) then

       call MAPL_AddExportSpec(GC,                                           &
          LONG_NAME      = 'surface_downwelling_shortwave_beam_flux_per_OBIO_band',&
          UNITS          = 'W m-2',                                          &
          SHORT_NAME     = 'DROBIO',                                         &
          DIMS           = MAPL_DimsHorzOnly,                                &
          UNGRIDDED_DIMS = (/ NB_OBIO /),                                    &
          VLOCATION      = MAPL_VLocationNone,                         __RC__)

       call MAPL_AddExportSpec(GC,                                           &
          LONG_NAME      = 'surface_downwelling_shortwave_diffuse_flux_per_OBIO_band',&
          UNITS          = 'W m-2',                                          &
          SHORT_NAME     = 'DFOBIO',                                         &
          DIMS           = MAPL_DimsHorzOnly,                                &
          UNGRIDDED_DIMS = (/ NB_OBIO /),                                    &
          VLOCATION      = MAPL_VLocationNone,                         __RC__)

    end if

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'cloud_area_fraction',                                   &
       UNITS      = '1' ,                                                    &
       SHORT_NAME = 'FCLD',                                                  &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                              __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'cloud_area_fraction_for_low_clouds',                    &
       UNITS      = '1' ,                                                    &
       SHORT_NAME = 'CLDLO',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'cloud_area_fraction_for_middle_clouds',                 &
       UNITS      = '1' ,                                                    &
       SHORT_NAME = 'CLDMD',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'cloud_area_fraction_for_high_clouds',                   &
       UNITS      = '1' ,                                                    &
       SHORT_NAME = 'CLDHI',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'total_cloud_area_fraction',                             &
       UNITS      = '1' ,                                                    &
       SHORT_NAME = 'CLDTT',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

#ifdef SOLAR_RADVAL
! Note: the four CLDxxSWHB diagnostics below represent super-layer cloud
! fractions based on essentially the same subcolumn cloud generation used
! by RRTMG SW but called from within the SOLAR UPDATE at the HEARTBEAT.
! They are GLOBAL (not just sunlit) fields and generated on the heartbeat.
! BUT, because subcolumn cloud generation is EXPENSIVE, asking for any of
! these exports will DOUBLE the cost of running the SOLAR GC. As such,
! they are for SPECIAL VALIDATION PURPOSES ONLY. No cost is incurred if
! they are not exported. But we encase them in SOLAR_RADVAL as an extra
! protection against their inadvertant use. Note, also, that they are NOT
! EXACTLY heartbeat versions of CLDxxSW, since they sample the heartbeat
! cloud fractions, not the less frequent snapshots used at REFRESH-frequency,
! and also since the generation inside UPDATE is on non-flipped vertical
! fields. This latter difference should be statistically insignificant.
! A re-coding to use vertically flipped fields as per RRTMG SW is possible
! but will be slightly slower, and was deemed unnecessary since the cloud
! fraction frequency difference will likely dominate.

    call MAPL_AddExportSpec(GC,                                              &
        SHORT_NAME = 'CLDTTSWHB',                                            &
        LONG_NAME  = 'total_cloud_area_fraction_rrtmg_sw_HEARTBEAT',         &
        UNITS      = '1',                                                    &
        DIMS       = MAPL_DimsHorzOnly,                                      &
        VLOCATION  = MAPL_VLocationNone,                              __RC__ )

    call MAPL_AddExportSpec(GC,                                              &
        SHORT_NAME = 'CLDHISWHB',                                            &
        LONG_NAME  = 'high-level_cloud_area_fraction_rrtmg_sw_HEARTBEAT',    &
        UNITS      = '1',                                                    &
        DIMS       = MAPL_DimsHorzOnly,                                      &
        VLOCATION  = MAPL_VLocationNone,                              __RC__ )

    call MAPL_AddExportSpec(GC,                                              &
        SHORT_NAME = 'CLDMDSWHB',                                            &
        LONG_NAME  = 'mid-level_cloud_area_fraction_rrtmg_sw_HEARTBEAT',     &
        UNITS      = '1',                                                    &
        DIMS       = MAPL_DimsHorzOnly,                                      &
        VLOCATION  = MAPL_VLocationNone,                              __RC__ )

    call MAPL_AddExportSpec(GC,                                              &
        SHORT_NAME = 'CLDLOSWHB',                                            &
        LONG_NAME  = 'low-level_cloud_area_fraction_rrtmg_sw_HEARTBEAT',     &
        UNITS      = '1',                                                    &
        DIMS       = MAPL_DimsHorzOnly,                                      &
        VLOCATION  = MAPL_VLocationNone,                              __RC__ )
#endif

    ! The TAUxx variants are ZERO when the super-layer is clear.
    ! These are the HISTORICAL exports, but are non-ideal as
    ! *in-cloud* diagnostics. The COTxx are to be prefered.
    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'in_cloud_optical_thickness_of_low_clouds',              &
       UNITS      = '1' ,                                                    &
       SHORT_NAME = 'TAULO',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'in_cloud_optical_thickness_of_middle_clouds',           &
       UNITS      = '1' ,                                                    &
       SHORT_NAME = 'TAUMD',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'in_cloud_optical_thickness_of_high_clouds',             &
       UNITS      = '1' ,                                                    &
       SHORT_NAME = 'TAUHI',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'in_cloud_optical_thickness_of_all_clouds__deprecated',  &
       UNITS      = '1' ,                                                    &
       SHORT_NAME = 'TAUTT',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'in_cloud_optical_thickness_of_all_clouds',              &
       UNITS      = '1' ,                                                    &
       SHORT_NAME = 'TAUTX',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    ! The COTxx variants are UNDEF when the super-layer is clear.
    ! They are preferred over TAUxx.
    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'in_cloud_optical_thickness_of_low_clouds_clrundef',     &
       UNITS      = '1' ,                                                    &
       SHORT_NAME = 'COTLO',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'in_cloud_optical_thickness_of_middle_clouds_clrundef',  &
       UNITS      = '1' ,                                                    &
       SHORT_NAME = 'COTMD',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'in_cloud_optical_thickness_of_high_clouds_clrundef',    &
       UNITS      = '1' ,                                                    &
       SHORT_NAME = 'COTHI',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'in_cloud_optical_thickness_of_all_clouds_clrundef',     &
       UNITS      = '1' ,                                                    &
       SHORT_NAME = 'COTTT',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    ! COT[DEN|NUM]xx allow a true cloud-fraction-weighted in-cloud optical
    ! thickness to be calculated via COTNUMxx / COTDENxx. Like the COTxx, clear
    ! values make no contribution, but unlike COTxx, each COT is weighted by a
    ! cloud fraction so that small clouds, which have a small radiative effect,
    ! get weighted accordingly. These provide the best estimate of radiatively
    ! effective in-cloud optical thicknesses, but require more advance post-
    ! processing (summing both NUM and DEN fields over the time-period required,
    ! and only then taking their quotient.)

    call MAPL_AddExportSpec(GC,                                                &
       LONG_NAME  = 'in_cloud_optical_thickness_of_low_clouds_denominator',    &
       UNITS      = '1' ,                                                      &
       SHORT_NAME = 'COTDENLO',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                         &
       VLOCATION  = MAPL_VLocationNone,                                  __RC__)

    call MAPL_AddExportSpec(GC,                                                &
       LONG_NAME  = 'in_cloud_optical_thickness_of_middle_clouds_denominator', &
       UNITS      = '1' ,                                                      &
       SHORT_NAME = 'COTDENMD',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                         &
       VLOCATION  = MAPL_VLocationNone,                                  __RC__)

    call MAPL_AddExportSpec(GC,                                                &
       LONG_NAME  = 'in_cloud_optical_thickness_of_high_clouds_denominator',   &
       UNITS      = '1' ,                                                      &
       SHORT_NAME = 'COTDENHI',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                         &
       VLOCATION  = MAPL_VLocationNone,                                  __RC__)

    call MAPL_AddExportSpec(GC,                                                &
       LONG_NAME  = 'in_cloud_optical_thickness_of_all_clouds_denominator',    &
       UNITS      = '1' ,                                                      &
       SHORT_NAME = 'COTDENTT',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                         &
       VLOCATION  = MAPL_VLocationNone,                                  __RC__)

    call MAPL_AddExportSpec(GC,                                                &
       LONG_NAME  = 'in_cloud_optical_thickness_of_low_clouds_numerator',      &
       UNITS      = '1' ,                                                      &
       SHORT_NAME = 'COTNUMLO',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                         &
       VLOCATION  = MAPL_VLocationNone,                                  __RC__)

    call MAPL_AddExportSpec(GC,                                                &
       LONG_NAME  = 'in_cloud_optical_thickness_of_middle_clouds_numerator',   &
       UNITS      = '1' ,                                                      &
       SHORT_NAME = 'COTNUMMD',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                         &
       VLOCATION  = MAPL_VLocationNone,                                  __RC__)

    call MAPL_AddExportSpec(GC,                                                &
       LONG_NAME  = 'in_cloud_optical_thickness_of_high_clouds_numerator',     &
       UNITS      = '1' ,                                                      &
       SHORT_NAME = 'COTNUMHI',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                         &
       VLOCATION  = MAPL_VLocationNone,                                  __RC__)

    call MAPL_AddExportSpec(GC,                                                &
       LONG_NAME  = 'in_cloud_optical_thickness_of_all_clouds_numerator',      &
       UNITS      = '1' ,                                                      &
       SHORT_NAME = 'COTNUMTT',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                         &
       VLOCATION  = MAPL_VLocationNone,                                  __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'in_cloud_optical_thickness_for_ice_clouds',             &
       UNITS      = '1' ,                                                    &
       SHORT_NAME = 'TAUCLI',                                                &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                              __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'in_cloud_optical_thickness_for_liquid_clouds',          &
       UNITS      = '1' ,                                                    &
       SHORT_NAME = 'TAUCLW',                                                &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                              __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'in_cloud_optical_thickness_for_falling_rain',           &
       UNITS      = '1' ,                                                    &
       SHORT_NAME = 'TAUCLR',                                                &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                              __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'in_cloud_optical_thickness_for_falling_snow',           &
       UNITS      = '1' ,                                                    &
       SHORT_NAME = 'TAUCLS',                                                &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                              __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'surface_net_downward_shortwave_flux_assuming_clear_sky',&
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'RSCS',                                                  &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'surface_net_downward_shortwave_flux',                   &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'RSRS',                                                  &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'surface_net_downward_shortwave_flux_assuming_clear_sky_and_no_aerosol',&
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'RSCSNA',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'surface_net_downward_shortwave_flux_assuming_no_aerosol',&
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'RSRSNA',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'surface_incoming_shortwave_flux',                       &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'SLRSF',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'surface_incoming_shortwave_flux_assuming_clear_sky',    &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'SLRSFC',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'surface_incoming_shortwave_flux_assuming_clean_sky',    &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'SLRSFNA',                                               &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'surface_incoming_shortwave_flux_assuming_clear_clean_sky',&
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'SLRSFCNA',                                              &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'surface_outgoing_shortwave_flux',                       &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'SLRSUF',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'surface_outgoing_shortwave_flux_assuming_clear_sky',    &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'SLRSUFC',                                               &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'surface_outgoing_shortwave_flux_assuming_clean_sky',    &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'SLRSUFNA',                                              &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'surface_outgoing_shortwave_flux_assuming_clear_clean_sky',&
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'SLRSUFCNA',                                             &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'toa_outgoing_shortwave_flux',                           &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'OSR',                                                   &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'toa_outgoing_shortwave_flux_assuming_clear_sky',        &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'OSRCLR',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'toa_outgoing_shortwave_flux_no_aerosol',                &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'OSRNA',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'toa_outgoing_shortwave_flux_no_aerosol__clear_sky',     &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'OSRCNA',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'toa_net_downward_shortwave_flux',                       &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'RSR',                                                   &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'toa_net_downward_shortwave_flux_assuming_clear_sky',    &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'RSC',                                                   &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'toa_net_downward_shortwave_flux_assuming_no_aerosol',   &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'RSRNA',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'toa_net_downward_shortwave_flux_assuming_clear_sky_and_no_aerosol',&
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'RSCNA',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'toa_incoming_shortwave_flux',                           &
       UNITS      = 'W m-2',                                                 &
       SHORT_NAME = 'SLRTP',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'surface_albedo',                                        &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'ALBEDO',                                                &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'surface_albedo_for_visible_beam',                       &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'ALBVR',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'surface_albedo_for_visible_diffuse',                    &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'ALBVF',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'surface_albedo_for_near_infrared_beam',                 &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'ALBNR',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'surface_albedo_for_near_infrared_diffuse',              &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'ALBNF',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    ! Three (one above now) different cos(SZA)s ...

    ! This one is instantaneous at the end of the UPDATE period,
    ! so it is consistent with the HISTORY files output time.
    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'cosine_of_the_solar_zenith_angle',                      &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'COSZ',                                                  &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    ! This one is the mean over the update period, so it will be half
    ! a heartbeat behind the HISTORY file time. But this is the flux-
    ! weighted and update-interval-averaged COSZ value used by the UPDATE.
    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'mean_cosine_of_the_solar_zenith_angle',                 &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'MCOSZ',                                                 &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                __RC__)

    call MAPL_AddExportSpec(GC,                                              &
        SHORT_NAME = 'CLDTMP',                                               &
        LONG_NAME  = 'cloud_top_temperature',                                &
        UNITS      = 'K',                                                    &
        DIMS       = MAPL_DimsHorzOnly,                                      &
        VLOCATION  = MAPL_VLocationNone,                               __RC__)

    call MAPL_AddExportSpec(GC,                                              &
        SHORT_NAME = 'CLDPRS',                                               &
        LONG_NAME  = 'cloud_top_pressure',                                   &
        UNITS      = 'Pa',                                                   &
        DIMS       = MAPL_DimsHorzOnly,                                      &
        VLOCATION  = MAPL_VLocationNone,                               __RC__)

!EOS

    ! Set Run method and use generic Initalize and Finalize methods
    call MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_RUN, Run, __RC__)
    call MAPL_GenericSetServices    (GC, __RC__)

    RETURN_(ESMF_SUCCESS)
  end subroutine SetServices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! !IROUTINE: RUN -- Run method for the SOLAR component

! !INTERFACE:

  subroutine RUN ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code
! /*
! !DESCRIPTION: Each time the Run method is called it fills all Exports for
!   which an allocated pointer is available. Exports are filled from the
!   normalized fluxes kept in the Internal state and the position of the Sun
!   for the current interval in the Clock. If MAPL's RunAlarm is ringing, it
!   also refreshes the normalized fluxes kept in the internal state by doing
!   a full transfer calculation valid for solar positions over a ``future
!   interval'' extending to the next anticipated ringing of the RunAlarm.
!   Whether this is done before or after the Exports are updated and the
!   exact definition of the ``future interval'' is controlled by a flag in
!   the configuration. \newline
!
!   A simple load balancing scheme is used that evens work between antipodal
!   processors. \newline
! \newline
! */

! !BUGS:
!
!\end{verbatim}
! \begin{itemize}
!  \item Deciding on the correct behavior for intermitent calls can be tricky.
!  \item Load-balancing communication needs to be upgraded to most up-to-date
!  ESMF machine model.
! \end{itemize}
!\begin{verbatim}
!

!EOP


! ErrLog Variables

    character(len=ESMF_MAXSTR) :: IAm
    character(len=ESMF_MAXSTR) :: COMP_NAME
    integer :: STATUS

! Local derived type aliases

    type (MAPL_MetaComp), pointer :: MAPL
    type (ESMF_Grid)              :: ESMFGRID
    type (ESMF_Config)            :: CF

    type (ty_RRTMGP_state), pointer :: rrtmgp_state => null()
    type (ty_RRTMGP_wrap)           :: wrap

! Local variables

    type (ESMF_Alarm)             :: ALARM
    type (ESMF_State)             :: INTERNAL
    type (ESMF_Time)              :: currentTime
    type (ESMF_TimeInterval)      :: intDT
    integer                       :: IM, JM, LM
    type (MAPL_SunOrbit)          :: ORBIT
    type (MAPL_VarSpec), pointer  :: ImportSpec(:)   => null()
    type (MAPL_VarSpec), pointer  :: ExportSpec(:)   => null()
    type (MAPL_VarSpec), pointer  :: InternalSpec(:) => null()
    real, pointer, dimension(:,:) :: LONS
    real, pointer, dimension(:,:) :: LATS

    real, pointer, dimension(:,:,:) :: ptr3d
    real, pointer, dimension(:,:  ) :: ptr2d

    type (ESMF_State)                     :: AERO
    character(len=ESMF_MAXSTR)            :: AS_FIELD_NAME
    integer                               :: AS_STATUS
    real, pointer,     dimension(:,:,:)   :: AS_PTR_3D
    real, pointer,     dimension(:,:,:)   :: AS_PTR_PLE
    real, pointer,     dimension(:,:,:)   :: AS_PTR_T
    real, pointer,     dimension(:,:,:)   :: AS_PTR_Q
    real, allocatable, dimension(:,:,:)   :: AS_ARR_RH
    real, allocatable, dimension(:,:,:)   :: AS_ARR_PL

    real, allocatable, dimension(:,:,:,:) :: AEROSOL_EXT
    real, allocatable, dimension(:,:,:,:) :: AEROSOL_SSA
    real, allocatable, dimension(:,:,:,:) :: AEROSOL_ASY

    integer :: band
    logical :: implements_aerosol_optics
    logical :: USE_RRTMGP, USE_RRTMGP_IRRAD
    logical :: USE_RRTMG,  USE_RRTMG_IRRAD
    logical :: USE_CHOU,   USE_CHOU_IRRAD
    integer :: NUM_BANDS_SOLAR, NUM_BANDS, TOTAL_RAD_BANDS
    logical :: SOLAR_TO_OBIO

    integer, parameter :: BANDS_SOLAR_OFFSET = 0

!PMN these three are redundant ... remove with zero-diff test later
    integer, parameter :: NB_CHOU  = 8         ! Num bands in SORAD calcs for Chou
    integer, parameter :: NB_RRTMG = 14        ! Num bands in SORAD calcs for RRTMG
    integer, parameter :: NB_RRTMGP = 14       ! Num bands in SORAD calcs for RRTMGP

    integer, parameter :: NB_CHOU_IRRAD  = 10  ! Num bands in IRRAD calcs for Chou
    integer, parameter :: NB_RRTMG_IRRAD = 16  ! Num bands in IRRAD calcs for RRTMG
    integer, parameter :: NB_RRTMGP_IRRAD = 16 ! Num bands in IRRAD calcs for RRTMGP

    integer :: CalledLast
    integer :: LCLDMH, LCLDLM
    integer :: YY, DOY
    integer :: K
    real    :: CO2
    real    :: PRS_LOW_MID
    real    :: PRS_MID_HIGH
    real    :: SC, HK(8), HK_IR_TEMP(3,10), HK_UV_TEMP(5), MG, SB
    integer :: SUNFLAG
    real, pointer, dimension(:) :: PREF

    logical :: REFRESH_FLUXES
    logical :: UPDATE_FIRST

    real, external                  :: getco2
    character(len=ESMF_MAXSTR)      :: MSGSTRING

    real, save :: CO2_0, SC_0, MG_0, SB_0
    data          CO2_0 /0.0/, SC_0/0.0/, MG_0/0.0/, SB_0/0.0/

    logical                    :: LoadBalance
    character(len=ESMF_MAXSTR) :: DYCORE
    integer                    :: SOLAR_LOAD_BALANCE
    integer                    :: SolarBalanceHandle
    integer                    :: MaxPasses

    character(len=ESMF_MAXPATHLEN) :: SolCycFileName
    logical                        :: USE_NRLSSI2
    logical                        :: PersistSolar

    logical :: do_no_aero_calc

    type(StringVector) :: string_vec
    type(StringVectorIterator) :: string_vec_iter
    character(len=:), pointer :: string_pointer

!=============================================================================

    ! Get the target components name and set-up traceback handle.
    call ESMF_GridCompGet (GC, name=COMP_NAME, GRID=ESMFGRID, __RC__ )
    Iam = trim(COMP_NAME) // "Run"

    ! Get my internal MAPL_Generic state
    call MAPL_GetObjectFromGC (GC, MAPL, __RC__)

    call MAPL_TimerOn (MAPL,"TOTAL"  ,__RC__)
    call MAPL_TimerOn (MAPL,"PRELIMS",__RC__)

    ! Get parameters from generic state.
    call MAPL_Get(MAPL,                                &
         IM                  = IM,                     &
         JM                  = JM,                     &
         LM                  = LM,                     &
         CF                  = CF,                     &
         LONS                = LONS,                   &
         LATS                = LATS,                   &
         RUNALARM            = ALARM,                  &
         ORBIT               = ORBIT,                  &
         INTERNALspec        = InternalSpec,           &
         IMPORTspec          = ImportSpec,             &
         EXPORTspec          = ExportSpec,             &
         INTERNAL_ESMF_STATE = INTERNAL,         __RC__)

    ! Get parameters from configuration
    call MAPL_GetResource (MAPL, PRS_LOW_MID,  'PRS_LOW_MID_CLOUDS:' , DEFAULT=70000., __RC__)
    call MAPL_GetResource (MAPL, PRS_MID_HIGH, 'PRS_MID_HIGH_CLOUDS:', DEFAULT=40000., __RC__)
    call MAPL_GetResource (MAPL, CO2,          'CO2:',                                 __RC__)
    call MAPL_GetResource (MAPL, SC,           'SOLAR_CONSTANT:',                      __RC__)
    call MAPL_GetResource (MAPL, SUNFLAG,      'SUN_FLAG:',            DEFAULT=0,      __RC__)

    ! Should we load balance solar radiation?
    ! For the single-column model, we always use the DATMO DYCORE.
    ! If this is our DYCORE, turn off load balancing.
    !---------------------------------------------
    call MAPL_GetResource (MAPL, DYCORE, 'DYCORE:', __RC__)
    call MAPL_GetResource (MAPL, SOLAR_LOAD_BALANCE, 'SOLAR_LOAD_BALANCE:', DEFAULT=1, __RC__)
    if (adjustl(DYCORE)=="DATMO" .OR. SOLAR_LOAD_BALANCE==0) then
       LoadBalance = .FALSE.
    else
       LoadBalance = .TRUE.
    end if

    ! Note: We set the default to 100 as that is the default in MAPL_LoadBalance which
    ! would have been used if not passed in
    call MAPL_GetResource (MAPL, MaxPasses, 'SOLAR_LB_MAX_PASSES:', DEFAULT=100, __RC__)

    ! Use time-varying co2
    call ESMF_ClockGet(CLOCK, currTIME=CURRENTTIME,       __RC__)
    call ESMF_TimeGet (CURRENTTIME, YY=YY, DayOfYear=DOY, __RC__)
    if(CO2<0.0) then
       CO2 = GETCO2(YY,DOY)
       write(MSGSTRING,'(A,I4,A,I3,A,e12.5)') &
            "Updated CO2 in solar for year/day ", YY, "/", DOY, " is ", CO2
       if( MAPL_AM_I_ROOT() ) then
          if( CO2_0.ne.CO2 ) then
             CO2_0  = CO2
             print *
             print *, trim(msgstring)
             print *
          endif
       endif
       call ESMF_LogWrite(MSGSTRING, ESMF_LOGMSG_INFO, __RC__)
    end if

    ! Decide which radiation to use:
    ! These USE_ flags are shared globally by contained SORADCORE() and Update_Flx()
    !-------------------------------------------------------------------------------
    call choose_solar_scheme (MAPL, &
      USE_RRTMGP,       USE_RRTMG,       USE_CHOU,       __RC__)
    call choose_irrad_scheme (MAPL, &
      USE_RRTMGP_IRRAD, USE_RRTMG_IRRAD, USE_CHOU_IRRAD, __RC__)

    ! Set number of solar bands
    if (USE_RRTMGP) then
      NUM_BANDS_SOLAR = NB_RRTMGP
    else if (USE_RRTMG) then
      NUM_BANDS_SOLAR = NB_RRTMG
    else
      NUM_BANDS_SOLAR = NB_CHOU
    end if

    ! Test to see if AGCM.rc is set up correctly for the Radiation selected
    !----------------------------------------------------------------------
    TOTAL_RAD_BANDS = NUM_BANDS_SOLAR
    if (USE_RRTMGP_IRRAD) then
      TOTAL_RAD_BANDS = TOTAL_RAD_BANDS + NB_RRTMGP_IRRAD
    else if (USE_RRTMG_IRRAD) then
      TOTAL_RAD_BANDS = TOTAL_RAD_BANDS + NB_RRTMG_IRRAD
    else
      TOTAL_RAD_BANDS = TOTAL_RAD_BANDS + NB_CHOU_IRRAD
    end if

    call MAPL_GetResource (MAPL, NUM_BANDS ,'NUM_BANDS:', __RC__)
    if (NUM_BANDS /= TOTAL_RAD_BANDS) then
      if (MAPL_am_I_Root()) then
         write (*,*) "NUM_BANDS is not set up correctly for the radiation combination selected:"
         write (*,*) "    SORAD RRTMG: ", USE_RRTMGP      , USE_RRTMG      , USE_CHOU
         write (*,*) "    IRRAD RRTMG: ", USE_RRTMGP_IRRAD, USE_RRTMG_IRRAD, USE_CHOU_IRRAD
         write (*,*) "Please check that your optics tables and NUM_BANDS are correct."
      end if
      _FAIL('Total number of radiation bands is inconsistent!')
   end if

   ! Decide if should make OBIO exports
   !-----------------------------------

    call MAPL_GetResource ( MAPL, DO_OBIO, Label="USE_OCEANOBIOGEOCHEM:",DEFAULT=0, RC=STATUS)
    VERIFY_(STATUS)

    SOLAR_TO_OBIO = (DO_OBIO/=0)

   ! Decide how to do solar forcing
   !-------------------------------

    call MAPL_GetResource (MAPL, SolCycFileName, "SOLAR_CYCLE_FILE_NAME:", DEFAULT='/dev/null', __RC__)
    if (SolCycFileName /= '/dev/null') THEN

       ! Solar forcing is from NRL SSI2 file for RRTMG[P].
       ! For chou-Suarez, the typical forcing is from internal tables, but a special
       ! file forcing is also possible.

       call MAPL_GetResource( MAPL, USE_NRLSSI2, "USE_NRLSSI2:", DEFAULT=.TRUE., __RC__)
       if (USE_NRLSSI2) then

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!! Use NRL SSI2 forcing file !!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         _ASSERT(USE_RRTMG .or. USE_RRTMGP, 'only RRTMG[P] can use NRLSSI2 currently')

         call MAPL_GetResource (MAPL, PersistSolar, "PERSIST_SOLAR:", DEFAULT=.TRUE., __RC__)
         call MAPL_SunGetSolarConstant (CLOCK, trim(SolCycFileName), &
            SC, MG, SB, PersistSolar=PersistSolar, __RC__)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! write(MSGSTRING,'(A,I4,A,I3,A,F8.3,A,F8.6,A,F9.4)') &                                          !
         !       "Solar Constants for year/day ", YY, "/", DOY, " are SC: ", SC, " MG: ", MG, " SB: ", SB !
         ! if( MAPL_AM_I_ROOT() ) then                                                                    !
         !    if( SC_0.ne.SC ) then                                                                       !
         !       SC_0  = SC                                                                               !
         !       MG_0  = MG                                                                               !
         !       SB_0  = SB                                                                               !
         !       print *                                                                                  !
         !       print *, trim(msgstring)                                                                 !
         !       print *                                                                                  !
         !    endif                                                                                       !
         ! endif                                                                                          !
         ! call ESMF_LogWrite (MSGSTRING, ESMF_LOGMSG_INFO, __RC__)                                       !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       else
         call MAPL_SunGetSolarConstant (CLOCK, trim(SolCycFileName), SC, HK=HK, __RC__)

         HK_UV_TEMP = HK(:5)

         do K=1,3
            HK_IR_TEMP(K,:)=HK_IR_OLD(K,:)*(HK(5+K)/sum(HK_IR_OLD(K,:)))
         end do
       end if

    else if (SC < 0.) then

       call MAPL_SunGetSolarConstant (CURRENTTIME, SC, HK, __RC__)

       HK_UV_TEMP = HK(:5)

       do K=1,3
          HK_IR_TEMP(K,:)=HK_IR_OLD(K,:)*(HK(5+K)/sum(HK_IR_OLD(K,:)))
       end do

       write(MSGSTRING,'(A,I4,A,I3,A,e12.5)') &
            "Updated Solar Constant for year/day ", YY, "/", DOY, " is ", SC
       if( MAPL_AM_I_ROOT() ) then
          if( SC_0.ne.SC ) then
             SC_0  = SC
             print *
             print *, trim(msgstring)
             print *
          endif
       endif
       call ESMF_LogWrite (MSGSTRING, ESMF_LOGMSG_INFO, __RC__)

    else
       HK_UV_TEMP = HK_UV_OLD
       HK_IR_TEMP = HK_IR_OLD
    end if

    ! Determine the model level separating high-middle and low-middle clouds
    !-----------------------------------------------------------------------

    ! Use the reference pressures to separate high, middle, and low clouds.
    call MAPL_GetPointer(IMPORT, PREF, 'PREF', __RC__)

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

    ! Determine calling sequence ...
    ! (CALLED_LAST == 1) => REFRESH called last, i.e., after UPDATE_EXPORT.
    ! This getresource is a kludge for now and needs to be fixed in the spec,
    ! because GC needs this info to know when to set the alarm, last or first
    ! step of interval. Right now it is always the last, which is only correct
    ! for called_last=1.
    !---------------------------

    call MAPL_GetResource (MAPL, CalledLast, 'CALLED_LAST:', default=1, __RC__)
    UPDATE_FIRST = CalledLast /= 0

    call MAPL_TimerOff(MAPL,"PRELIMS",__RC__)

    ! Update the Sun position and weight the export variables
    ! -------------------------------------------------------
    if (UPDATE_FIRST) then
       call MAPL_TimerOn  (MAPL,"UPDATE",__RC__)
       call UPDATE_EXPORT (IM,JM,LM,__RC__)
       call MAPL_TimerOff (MAPL,"UPDATE",__RC__)
    end if

    ! Periodically, refresh the internal state with a full solar calc
    ! ---------------------------------------------------------------
    REFRESH_FLUXES = ESMF_AlarmIsRinging (ALARM, __RC__)

    REFRESH: if (REFRESH_FLUXES) then
       call MAPL_TimerOn (MAPL,"REFRESH",__RC__)

       call ESMF_AlarmRingerOff (ALARM, __RC__)
       call ESMF_ClockGet (CLOCK, currTIME=CURRENTTIME, __RC__)

       ! Beginning of REFRESH period is current time PLUS offset intDT
       ! -------------------------------------------------------------
       if (UPDATE_FIRST) then
          ! The UPDATE is already done, so the REFRESH interval should start one
          ! timestep beyond current time so it is consistent with the NEXT update.
          call ESMF_ClockGet(CLOCK, timeSTEP=intDT, __RC__)
       else
          ! The UPDATE will occur after the REFRESH, so both update and refresh
          ! periods should begin at the current time.
          call ESMF_TimeIntervalSet(intDT, s=0, __RC__)
       end if

       ! Get optical properties of radiatively active aerosols
       ! -----------------------------------------------------
       call MAPL_TimerOn(MAPL,"-AEROSOLS",__RC__)
       call ESMF_StateGet(IMPORT,'AERO',AERO,__RC__)
       call ESMF_AttributeGet(AERO, &
          name='implements_aerosol_optics_method', &
          value=implements_aerosol_optics,__RC__)
       if (implements_aerosol_optics) then

           ! set RH for aerosol optics
           call ESMF_AttributeGet(AERO, &
              name='relative_humidity_for_aerosol_optics', &
              value=AS_FIELD_NAME,__RC__)
           if (AS_FIELD_NAME /= '') then
              call MAPL_GetPointer(IMPORT,AS_PTR_PLE,'PLE',__RC__)
              call MAPL_GetPointer(IMPORT,AS_PTR_Q,  'QV', __RC__)
              call MAPL_GetPointer(IMPORT,AS_PTR_T,  'T',  __RC__)
              allocate(AS_ARR_RH(IM,JM,LM),AS_ARR_PL(IM,JM,LM),__STAT__)
              AS_ARR_PL = 0.5 * (AS_PTR_PLE(:,:,1:LM) + AS_PTR_PLE(:,:,0:LM-1))
              AS_ARR_RH = AS_PTR_Q / MAPL_EQSAT(AS_PTR_T,PL=AS_ARR_PL)
              call MAPL_GetPointer(AERO,AS_PTR_3D,trim(AS_FIELD_NAME),__RC__)
              AS_PTR_3D = AS_ARR_RH
              deallocate(AS_ARR_RH,AS_ARR_PL,__STAT__)
           end if

           ! set PLE for aerosol optics
           call ESMF_AttributeGet(AERO, &
              name='air_pressure_for_aerosol_optics', &
              value=AS_FIELD_NAME,__RC__)
           if (AS_FIELD_NAME /= '') then
              call MAPL_GetPointer(IMPORT,AS_PTR_PLE,'PLE',__RC__)
              call MAPL_GetPointer(AERO,AS_PTR_3D,trim(AS_FIELD_NAME),__RC__)
              AS_PTR_3D = AS_PTR_PLE
           end if

           ! allocate memory for TOTAL aerosol ext, ssa and asy at all solar bands
           allocate(AEROSOL_EXT(IM,JM,LM,NUM_BANDS_SOLAR), &
                    AEROSOL_SSA(IM,JM,LM,NUM_BANDS_SOLAR), &
                    AEROSOL_ASY(IM,JM,LM,NUM_BANDS_SOLAR), __STAT__)

           ! zero by default
           ! (in case aero provider cant provide some of them)
           AEROSOL_EXT = 0.
           AEROSOL_SSA = 0.
           AEROSOL_ASY = 0.

           ! compute aerosol optics at all solar bands
           SOLAR_BANDS: do band = 1, NUM_BANDS_SOLAR
              call ESMF_AttributeSet(AERO, &
                 name='band_for_aerosol_optics', &
                 value=(BANDS_SOLAR_OFFSET+band),__RC__)

              ! execute the aero provider's optics method
              call ESMF_MethodExecute(AERO, &
                 label="run_aerosol_optics", &
                 userRC=AS_STATUS, RC=STATUS)
              VERIFY_(AS_STATUS)
              VERIFY_(STATUS)

              ! EXT from AERO_PROVIDER
              call ESMF_AttributeGet(AERO, &
                 name='extinction_in_air_due_to_ambient_aerosol', &
                 value=AS_FIELD_NAME,__RC__)
              if (AS_FIELD_NAME /= '') then
                 call MAPL_GetPointer(AERO,AS_PTR_3D,trim(AS_FIELD_NAME),__RC__)
                 if (associated(AS_PTR_3D)) AEROSOL_EXT(:,:,:,band) = AS_PTR_3D
              end if

              ! SSA from AERO_PROVIDER (actually EXT * SSA)
              call ESMF_AttributeGet(AERO, &
                 name='single_scattering_albedo_of_ambient_aerosol', &
                 value=AS_FIELD_NAME,__RC__)
              if (AS_FIELD_NAME /= '') then
                 call MAPL_GetPointer(AERO,AS_PTR_3D,trim(AS_FIELD_NAME),__RC__)
                 if (associated(AS_PTR_3D)) AEROSOL_SSA(:,:,:,band) = AS_PTR_3D
              end if

              ! ASY from AERO_PROVIDER (actually EXT * SSA * ASY)
              call ESMF_AttributeGet(AERO, &
                 name='asymmetry_parameter_of_ambient_aerosol', &
                 value=AS_FIELD_NAME,__RC__)
              if (AS_FIELD_NAME /= '') then
                 call MAPL_GetPointer(AERO,AS_PTR_3D,trim(AS_FIELD_NAME),__RC__)
                 if (associated(AS_PTR_3D)) AEROSOL_ASY(:,:,:,band) = AS_PTR_3D
              end if

           end do SOLAR_BANDS

       end if  ! implements_aerosol_optics
       call MAPL_TimerOff(MAPL,"-AEROSOLS", __RC__)

       ! Optional without-aerosol diagnostics
       ! ------------------------------------

       ! are without-aerosol exports requested?
       do_no_aero_calc = .false.

       call string_vec%push_back('FSWNA')
       call string_vec%push_back('FSWUNA')
       call string_vec%push_back('FSWDNA')
       call string_vec%push_back('FSCNA')
       call string_vec%push_back('FSCUNA')
       call string_vec%push_back('FSCDNA')
       call string_vec%push_back('FSWBANDNA')

       string_vec_iter = string_vec%begin()
       do while ( string_vec_iter /= string_vec%end() )
          string_pointer => string_vec_iter%get()
          call MAPL_GetPointer( EXPORT, ptr3d, string_pointer, __RC__)
          do_no_aero_calc = (do_no_aero_calc .or. associated(ptr3d))
          call string_vec_iter%next()
       end do

       if (.not. do_no_aero_calc) then

         call string_vec%clear()

         call string_vec%push_back('RSRNA')
         call string_vec%push_back('RSRSNA')
         call string_vec%push_back('OSRNA')
         call string_vec%push_back('RSCNA')
         call string_vec%push_back('RSCSNA')
         call string_vec%push_back('OSRCNA')
         call string_vec%push_back('SLRSFNA')
         call string_vec%push_back('SLRSUFNA')
         call string_vec%push_back('SLRSFCNA')
         call string_vec%push_back('SLRSUFCNA')

         string_vec_iter = string_vec%begin()

         do while ( string_vec_iter /= string_vec%end() )
            string_pointer => string_vec_iter%get()
            call MAPL_GetPointer( EXPORT, ptr2d, string_pointer, __RC__)
            do_no_aero_calc = (do_no_aero_calc .or. associated(ptr2d))
            call string_vec_iter%next()
         end do
       end if

       if (do_no_aero_calc) then

          ! do a calculation without aerosols:
          !   this just sets the no-aerosol internals ---
          !   the exports are derived from the internals in update_export()
          call SORADCORE(IM,JM,LM,           &
               include_aerosols = .false.,   &
               CURRTIME = CURRENTTIME+intDT, &
               MaxPasses = MaxPasses,        &
               LoadBalance = LoadBalance,    &
               __RC__)
       else

          ! otherwise, zero the no-aerosol internals
          call string_vec%clear()

          call string_vec%push_back('FSWNAN')
          call string_vec%push_back('FSWUNAN')
          call string_vec%push_back('FSCNAN')
          call string_vec%push_back('FSCUNAN')
          call string_vec%push_back('FSWBANDNAN')
          string_vec_iter = string_vec%begin()
          do while ( string_vec_iter /= string_vec%end() )
             string_pointer => string_vec_iter%get()
             call MAPL_GetPointer( INTERNAL, ptr3d, string_pointer, __RC__)
             ptr3d = 0.
             call string_vec_iter%next()
          end do

       end if

       ! Regular with-aerosol calculations
       ! ---------------------------------
       call SORADCORE(IM,JM,LM,                    &
                      include_aerosols = .true.,   &
                      CURRTIME = CURRENTTIME+intDT,&
                      MaxPasses = MaxPasses,       &
                      LoadBalance = LoadBalance,   &
                      __RC__)

       ! Clean up aerosol optical properties
       ! -----------------------------------
       if (implements_aerosol_optics) then
          deallocate(AEROSOL_EXT,__STAT__)
          deallocate(AEROSOL_SSA,__STAT__)
          deallocate(AEROSOL_ASY,__STAT__)
       end if

       call MAPL_TimerOff(MAPL,"REFRESH",__RC__)
    endif REFRESH

    ! Update the Sun position and weight the export variables
    ! -------------------------------------------------------
    if (.not.UPDATE_FIRST) then
       call MAPL_TimerOn  (MAPL,"UPDATE",__RC__)
       call UPDATE_EXPORT (IM,JM,LM,     __RC__)
       call MAPL_TimerOff (MAPL,"UPDATE",__RC__)
    end if

    call MAPL_TimerOff (MAPL,"TOTAL",__RC__)
    RETURN_(ESMF_SUCCESS)

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine SORADCORE(IM,JM,LM,include_aerosols,CURRTIME,MaxPasses,LoadBalance,RC)

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
      use mo_fluxes_byband,           only: ty_fluxes_byband
      use mo_rte_sw,                  only: rte_sw
      use mo_load_coefficients,       only: load_and_init
      use mo_load_cloud_coefficients, only: load_cld_lutcoeff, load_cld_padecoeff

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
#endif

      ! for RRTMGP (use implicit inside RRTMG)
      use cloud_condensate_inhomogeneity, only: condensate_inhomogeneous, zcw_lookup
      use cloud_subcol_gen, only : &
        correlation_length_cloud_fraction, correlation_length_condensate

      implicit none

      integer,           intent(IN ) :: IM, JM, LM
      logical,           intent(IN ) :: include_aerosols
      type (ESMF_Time),  intent(IN ) :: CURRTIME
      integer,           intent(IN ) :: MaxPasses
      logical,           intent(IN ) :: LoadBalance
      integer, optional, intent(OUT) :: RC

      ! Locals

      character(len=ESMF_MAXSTR)     :: IAm
      integer                        :: STATUS

      type (ESMF_TimeInterval)       :: TINT
      type (ESMF_DELayout)           :: LAYOUT
      type (ESMF_Array)              :: ARRAY
      type (ESMF_FieldBundle)        :: BUNDLE
      type (ESMF_Field)              :: FIELD

      type (ESMF_VM)                 :: VM
      integer                        :: COMM

      real,    dimension(IM,JM)      :: ZTH, SLR
      logical, dimension(IM,JM)      :: daytime

      ! Daytime ONLY copy of variables
      ! ------------------------------

      ! inputs
      real, pointer, dimension(:,:)  :: PLE, CH4, N2O, T, Q, OX, CL, &
                                        QL, QI, QR, QS, RL, RI, RR, RS
      real, pointer, dimension(:)    :: TS, ALBNR, ALBNF, ALBVR, ALBVF, &
                                        Ig1D, Jg1D, ALAT, SLR1D, ZT

      ! outputs (via internals)
      real, pointer, dimension(:,:)  :: FSW, FSC, FSWA, FSCA, FSWU, FSCU, FSWUA, FSCUA, &
                                        FSWBAND, FSWBANDA
      real, pointer, dimension(:)    :: UVRR, UVRF, PARR, PARF, NIRR, NIRF

      ! outputs used for OBIO
      real, pointer, dimension(:,:)  :: DRBAND, DFBAND

      ! REFRESH exports (via internals)
      real, pointer, dimension(:)    :: COSZSW
      real, pointer, dimension(:)    :: CLDTS, CLDHS, CLDMS, CLDLS, &
#ifdef SOLAR_RADVAL
                                        TAUTP, TAUHP, TAUMP, TAULP, &
#endif
                                        COTTP, COTHP, COTMP, COTLP, &
                                        COTDTP, COTDHP, COTDMP, COTDLP, &
                                        COTNTP, COTNHP, COTNMP, COTNLP

#ifdef SOLAR_RADVAL
      real, pointer, dimension(:)    :: CDSDTP, CDSDHP, CDSDMP, CDSDLP, &
                                        CDSNTP, CDSNHP, CDSNMP, CDSNLP, &

                                        COTLDTP, COTLDHP, COTLDMP, COTLDLP, &
                                        COTLNTP, COTLNHP, COTLNMP, COTLNLP, &
                                        COTIDTP, COTIDHP, COTIDMP, COTIDLP, &
                                        COTINTP, COTINHP, COTINMP, COTINLP, &

                                        CDSLDTP, CDSLDHP, CDSLDMP, CDSLDLP, &
                                        CDSLNTP, CDSLNHP, CDSLNMP, CDSLNLP, &
                                        CDSIDTP, CDSIDHP, CDSIDMP, CDSIDLP, &
                                        CDSINTP, CDSINHP, CDSINMP, CDSINLP, &

                                        SSALDLP, SSALNLP, SSAIDLP, SSAINLP, &
                                        SSALDMP, SSALNMP, SSAIDMP, SSAINMP, &
                                        SSALDHP, SSALNHP, SSAIDHP, SSAINHP, &
                                        SSALDTP, SSALNTP, SSAIDTP, SSAINTP, &

                                        SDSLDLP, SDSLNLP, SDSIDLP, SDSINLP, &
                                        SDSLDMP, SDSLNMP, SDSIDMP, SDSINMP, &
                                        SDSLDHP, SDSLNHP, SDSIDHP, SDSINHP, &
                                        SDSLDTP, SDSLNTP, SDSIDTP, SDSINTP, &

                                        ASMLDLP, ASMLNLP, ASMIDLP, ASMINLP, &
                                        ASMLDMP, ASMLNMP, ASMIDMP, ASMINMP, &
                                        ASMLDHP, ASMLNHP, ASMIDHP, ASMINHP, &
                                        ASMLDTP, ASMLNTP, ASMIDTP, ASMINTP, &

                                        ADSLDLP, ADSLNLP, ADSIDLP, ADSINLP, &
                                        ADSLDMP, ADSLNMP, ADSIDMP, ADSINMP, &
                                        ADSLDHP, ADSLNHP, ADSIDHP, ADSINHP, &
                                        ADSLDTP, ADSLNTP, ADSIDTP, ADSINTP, &

                                        FORLDLP, FORLNLP, FORIDLP, FORINLP, &
                                        FORLDMP, FORLNMP, FORIDMP, FORINMP, &
                                        FORLDHP, FORLNHP, FORIDHP, FORINHP, &
                                        FORLDTP, FORLNTP, FORIDTP, FORINTP
#endif

      ! variables for RRTMG code
      ! ------------------------

      integer :: ICEFLGSW        ! Flag for ice particle specification
      integer :: LIQFLGSW        ! Flag for liquid droplet specification

      real,    allocatable, dimension(:,:)   :: TLEV, TLEV_R, PLE_R
      real,    allocatable, dimension(:,:)   :: FCLD_R, CLIQWP, CICEWP, RELIQ, REICE
      real,    allocatable, dimension(:,:,:) :: TAUAER, SSAAER, ASMAER
      real,    allocatable, dimension(:,:)   :: DPR, PL_R, ZL_R, T_R, Q_R, O2_R, O3_R, CO2_R, CH4_R

      integer, allocatable, dimension(:,:)   :: CLEARCOUNTS
      real,    allocatable, dimension(:,:)   :: SWUFLX,  SWDFLX,  SWUFLXC,  SWDFLXC
      real,    allocatable, dimension(:,:)   :: SWUFLXR, SWDFLXR, SWUFLXCR, SWDFLXCR

      ! pmn: should we update these?
      real, parameter :: O2   = 0.2090029E+00 ! preexisting
      real, parameter :: N2   = 0.7906400E+00 ! approx from rrtmgp input file
      real, parameter :: CO   = 0.0           ! currently zero

      real    :: ADJES, DIST
      integer :: DYOFYR
      integer :: NCOL
      integer :: RPART, IAER, NORMFLX

      integer                   :: ISOLVAR
      real, dimension(2)        :: INDSOLVAR
      real, dimension(nb_rrtmg) :: BNDSOLVAR
      real                      :: SOLCYCFRAC

      ! variables for RRTMGP code
      ! -------------------------

      ! conversion factor (see below)
      real(wp), parameter :: cwp_fac = real(1000./MAPL_GRAV,kind=wp)

      ! gpoint limits for each band (2,nbnd)
      integer, dimension(:,:),      allocatable         :: band_lims_gpt

      ! solar inputs: (ncol) and (nbnd,ncol)
      real(wp), dimension(:),       allocatable         :: tsi, mu0
      real(wp), dimension(:,:),     allocatable         :: sfc_alb_dir, sfc_alb_dif

      ! per g-point toa flux (ncols_block,ngpt) [W/m2 NORMAL to solar beam]
      real(wp), dimension(:,:),     allocatable         :: toa_flux

      ! input arrays: dimensions (ncol,nlay[+1]) [Pa,K]
      real(wp), dimension(:,:),     allocatable         :: dummy_wp
      real(wp), dimension(:,:),     allocatable         :: p_lay, t_lay, dp_wp
      real(wp), dimension(:,:),     allocatable         :: p_lev

      ! inter-layer separations (from mid-points) (ncol,nlay-1) [m]
      real(wp), dimension(:,:),     allocatable         :: dzmid

      ! fluxes that we actually need
      ! NB: fluxes_byand makes available fluxes%[bnd_]flux_[up|dn|net|dn_dir->"dir"].
      real(wp), dimension(:),       allocatable         :: flux_dn_top
      real(wp), dimension(:,:),     allocatable, target :: flux_up_clrsky, flux_net_clrsky
      real(wp), dimension(:,:),     allocatable, target :: flux_up_allsky, flux_net_allsky
      real(wp), dimension(:,:,:),   allocatable, target :: bnd_flux_dn_allsky, bnd_flux_net_allsky, &
                                                           bnd_flux_dir_allsky

      ! derived types for interacting with RRTMGP
      ! (cloud_optics generates cloud_props_bnd from loaded
      ! coefficients and cloud physical properties)
      type(ty_gas_optics_rrtmgp), pointer               :: k_dist
      type(ty_gas_concs)                                :: gas_concs, gas_concs_block
      type(ty_cloud_optics_rrtmgp)                      :: cloud_optics
      type(ty_fluxes_byband)                            :: fluxes_clrsky, fluxes_allsky

      ! PMN: my earlier RRTMGP implementations used cloud_props for liq and ice combined,
      ! but now, to allow separate delta-scaling for the two phases, we keep separate liq and
      ! ice properties, and combine them later. There may be some speedup possible here, but 
      ! to allow for future more independent phases (e.g., separate condensate inhomogeneity
      ! for the phases), we keep the phase optical properties separate as long as possible.

      ! The band-space (ncol,nlay,nbnd) aerosol and in-cloud optical properties
      ! Polymorphic with dynamic type (#streams) defined later
      class(ty_optical_props_arry), allocatable :: aer_props
      class(ty_optical_props_arry), allocatable :: cloud_props_bnd_liq, cloud_props_bnd_ice

      ! The g-point cloud optical properties used for mcICA
      class(ty_optical_props_arry), allocatable :: cloud_props_gpt_liq, cloud_props_gpt_ice

      ! The g-point optical properties used in RT calculations
      ! Polymorphic with dynamic type (#streams) defined later
      class(ty_optical_props_arry), allocatable :: optical_props

      ! RRTMGP locals
      logical :: top_at_1, partial_block, need_aer_optical_props
      logical :: gen_mro, cond_inhomo
      logical :: rrtmgp_delta_scale, rrtmgp_use_rrtmg_iceflg3_like_forwice
      integer :: nbnd, ngpt, nmom, icergh
      integer :: ib, b, nBlocks, colS, colE, ncols_block, &
                 partial_blockSize, icol, isub, ilay, igpt
      real(wp), allocatable :: t_lev(:) ! (ncol)
      character(len=ESMF_MAXPATHLEN) :: k_dist_file, cloud_optics_file
      character(len=ESMF_MAXSTR)     :: error_msg
      character(len=128)             :: cloud_optics_type, cloud_overlap_type
      type (ESMF_Time)               :: ReferenceTime
      type (ESMF_TimeInterval)       :: RefreshInterval
      real :: cld_frac, sigma_qcw, wgt
      real :: stautp, stauhp, staump, staulp
      real :: sltautp, sltauhp, sltaump, sltaulp
      real :: sitautp, sitauhp, sitaump, sitaulp
#ifdef SOLAR_RADVAL
      real :: sltaussatp, sltaussahp, sltaussamp, sltaussalp
      real :: sitaussatp, sitaussahp, sitaussamp, sitaussalp
      real :: sltaussagtp, sltaussaghp, sltaussagmp, sltaussaglp
      real :: sitaussagtp, sitaussaghp, sitaussagmp, sitaussaglp
      real :: sltaussaftp, sltaussafhp, sltaussafmp, sltaussaflp
      real :: sitaussaftp, sitaussafhp, sitaussafmp, sitaussaflp
#endif

      ! radice interpolation for forwice
      integer :: radidx
      real(wp) :: radice_lwr, radice_upr, radice
      real(wp) :: radfac, rfint, fdelta

      ! for global gcolumn index seeding of PRNGs
      integer :: iBeg, iEnd, jBeg, jEnd
      integer :: IM_World, JM_World, Gdims(3)
      integer, dimension(IM,JM) :: Ig, Jg

      ! a column random number generator
#ifdef HAVE_MKL
      type(ty_rng_mklvsl_plus) :: rng
#endif
      integer, dimension(:), allocatable :: seeds

      ! uniform random numbers needed by mcICA (ngpt,nlay,rrtmgp_blocksize)
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

      ! forward scattering fraction for cloud droplets and ice crystals (ncols_block,nlay,ngpt)
      real(wp), dimension(:,:,:), allocatable :: forwliq, forwice

      ! TEMP ... see below
      real(wp) :: press_ref_min, ptop
      real(wp) ::  temp_ref_min, tmin
      real(wp), parameter :: ptop_increase_OK_fraction = 0.01_wp
      real(wp) :: tmin_increase_OK_Kelvin

      ! block size for efficient column processing (set from resource file)
      integer :: rrtmgp_blockSize

      ! For Aerosol
      real, pointer, dimension(:,:,:)     :: taua, ssaa, asya
      real, pointer, dimension(:,:,:)     :: BUFIMP_AEROSOL_EXT => null()
      real, pointer, dimension(:,:,:)     :: BUFIMP_AEROSOL_SSA => null()
      real, pointer, dimension(:,:,:)     :: BUFIMP_AEROSOL_ASY => null()
      real, allocatable, dimension(:,:,:) :: BUF_AEROSOL

      ! LoadBalance and general
      integer :: I, J, K, L, i1, iN
      integer, pointer :: pi1, piN
      integer, target :: i1Out, iNOut, i1InOut, iNInOut
      real, pointer :: QQ3(:,:,:), RR3(:,:,:), ptr3(:,:,:), ptr4(:,:,:,:)
      real, pointer :: ptr2(:,:), RH(:,:), PL(:,:), O3(:,:), PLhPa(:,:)
      integer :: dims, NumLit, Num2do, num_aero_vars
      character(len=ESMF_MAXSTR) :: short_name
      integer, pointer :: ugdims(:) => null()
      logical, allocatable :: IntInOut(:)
      character(len=ESMF_MAXSTR), allocatable :: NamesInp(:), NamesInt(:)
      integer, allocatable :: SlicesInp(:), SlicesInt(:)
      real, target, allocatable :: BufInp(:), BufInOut(:), BufOut(:)
      integer, allocatable :: rgDim(:), ugDim(:)
      real, pointer :: buf(:)
      integer :: NumImp, NumInt, NumInp
      integer :: NumMax, HorzDims(2)
      integer :: ibinary
      real :: def

      IAm = trim(COMP_NAME)//"Soradcore"
      call MAPL_TimerOn(MAPL,"-MISC")

! Get the average insolation for the next alarm "REFRESH" interval
!-----------------------------------------------------------------
! @ In standard (legacy) mode, this longer REFRESH interval forms the basis of
! a normalized full solar calculation that is simply scaled at each hearbeat (in
! UPDATE_EXPORTS) by the TOA projected solar input. This scaled update is extremely
! quick, but it lacks some important aspects of the full calculation, namely the
! pathlength (cf. projection) effect of the updated solar position, and all the
! changes caused by variations in surface albedo and atmospheric properties within
! the REFRESH period.

      call ESMF_AlarmGet(ALARM, RINGINTERVAL=TINT, __RC__)
      call MAPL_SunGetInsolation(  &
              LONS, LATS,          &
              ORBIT, ZTH, SLR,     &
              INTV = TINT,         &
              currTime = currTime, &
              TIME = SUNFLAG,      &
              DIST = DIST,         &
              __RC__)

      ! convert SLR to an ABSOLUTE downward flux [W/m2] for normalization purposes later
      ! (SLR from MAPL_SunGetInsolation() already contains the DIST and ZTH effects)
      SLR = SLR * SC

      ! prepare global gridcolumn indicies needed by random number generators
      ! get indicies of local rectangular grid
      call MAPL_GridGet(ESMFGRID, globalCellCountPerDim=Gdims, __RC__)
      IM_World = Gdims(1); JM_World = Gdims(2)
      call MAPL_GridGetInterior (ESMFGRID,iBeg,iEnd,jBeg,jEnd)
      do J=1,JM
        do I=1,IM
          Ig(I,J) = iBeg + I - 1
          Jg(I,J) = jBeg + J - 1
        end do
      end do

      call MAPL_TimerOff(MAPL,"-MISC")

!  Load balancing by packing the lit points and sharing work with night regions
!------------------------------------------------------------------------------

      call MAPL_TimerOn(MAPL,"-BALANCE")

!  Identify lit soundings with the daytime mask
!----------------------------------------------

!  The load balancer does not work if there are no lit points. This is only
!  important model-wise with the single-column model. Note we must protect
!  ZTH since in solar, we divide by ZTH and, thus, we will get a divide-by-
!  zero if not protected.
!--------------------------------------------------------------------------

      if (adjustl(DYCORE)=="DATMO") ZTH = max(.0001,ZTH)

      daytime = ZTH > 0.
      NumLit  = count(daytime)

!  Create a balancing strategy. This is a collective call on the communicator
!  of the current VM. The original, unbalanced local work consists of (OrgLen)
!  NumLit soundings, which may be zero. The local work after implementing the
!  strategy consists of (BalLen) Num2do soundings, which is generally non-zero.
!  The data movement to implement this strategy will occur when MAPL_BalanceWork
!  is called to "distribute" excess work to less busy processors and later to
!  "retrieve" that work to its home processor. Because the data balancing will be
!  done "in place", the in-out buffer must be large enough to accomodate the data
!  held at each stage (pass) of the balancing. This can be larger than the max of
!  the initial and final sizes; so the required size is passed back in BufLen.
!------------------------------------------------------------------------------------

      call ESMF_VMGetCurrent(VM, __RC__)
      call ESMF_VMGet(VM, mpiCommunicator=COMM, __RC__)

      call MAPL_TimerOn(MAPL,"--CREATE")

      if (LoadBalance) then
         call MAPL_BalanceCreate( &
            OrgLen=NumLit, Comm=COMM, MaxPasses=MaxPasses, Handle=SolarBalanceHandle, &
            BalLen=Num2do, BufLen=NumMax, __RC__)
      else
         Num2do = NumLit
         NumMax = NumLit
      end if

      call MAPL_TimerOff(MAPL,"--CREATE")

!  The number of Input and Output/InOut variables to the load balancing.
!    The Input number is five more than the number of IMPORTS because the
!  component needs the LATS, SLR and ZTH from MAPL and the global grid-
!  column indicies Ig and Jg.
!    The Outputs and InOuts are all INTERNAL variables.
!--------------------------------------------------------------

      NumImp = size(ImportSpec)
      NumInt = size(InternalSpec)

      ! Inputs to load balancing:
      ! All imports plus Ig, Jg, LATS, SLR & ZTH.
      ! Vertical only imports (PREF) are explicitly skipped later.
      NumInp = NumImp + 5

      allocate( &
         SlicesInp(NumInp), NamesInp(NumInp), &
         SlicesInt(NumInt), NamesInt(NumInt), &
         IntInOut(NumInt), rgDim(NumInt), ugDim(NumInt), &
         __STAT__)

      HorzDims = (/IM,JM/)

      ! num_aero_vars for aerosol optics from aerosol bundle
      if (include_aerosols .and. implements_aerosol_optics) then
         num_aero_vars = 3
      else
         num_aero_vars = 0
      end if

! @@@@@@@@@@@@@@
! @@@ Inputs @@@
! @@@@@@@@@@@@@@

      ! Count the 2D slices in the Input variables and calculate
      ! the required length of the 1D Input buffer (BufInp)
      ! --------------------------------------------------------
      INPUT_VARS_1: do k=1,NumInp

         ! Get names and dimensions of Inputs
         if (k <= NumImp) then
            call MAPL_VarSpecGet(ImportSpec(k), &
               DIMS=dims, SHORT_NAME=NamesInp(k), __RC__)
         else
            dims = MAPL_DIMSHORZONLY
            if      (k == NumImp+1) then
               NamesInp(k) = "Ig"
            else if (k == NumImp+2) then
               NamesInp(k) = "Jg"
            else if (k == NumImp+3) then
               NamesInp(k) = "LATS"
            else if (k == NumImp+4) then
               NamesInp(k) = "SLR"
            else if (k == NumImp+5) then
               NamesInp(k) = "ZTH"
            end if
         end if

         ! Skip vertical only inputs (PREF). They dont require
         ! load-balancing since they have no horizontal dimension.
         if (dims == MAPL_DIMSVERTONLY) then
            SlicesInp(k) = 0
            cycle
         end if

         ! If Import is aerosol bundle, make space for aerosol optical props.
         ! Note: assume all aerosol species are dimensioned by LM levels.
         !    This will be asserted later.

         if (NamesInp(k) == "AERO") then

            SlicesInp(k) = LM * num_aero_vars * NUM_BANDS_SOLAR

         else  ! Non-aerosol input

            select case(dims)
               case(MAPL_DIMSHORZVERT)
                  ! We currently assume this case is 3D
                  call ESMFL_StateGetPointerToData(IMPORT,ptr3,NamesInp(k),__RC__)
                  SlicesInp(k) = size(ptr3,3)

               case(MAPL_DIMSHORZONLY)
                  SlicesInp(k) = 1

               case default
                  _FAIL('invalid dimension for SOLAR import')
            end select

         end if

      enddo INPUT_VARS_1

      ! Allocate buffer with enough space to hold both the unbalanced
      ! and balanced data on the local PE for Input vars. The inner
      ! dimension of its 2D representation must be NumMax.
      ! -------------------------------------------------------------
      allocate(BufInp(NumMax*sum(SlicesInp)),__STAT__)
      BufInp = MAPL_UNDEF

      ! Loop over imports, packing into the buffer that will be
      ! load balanced and used in the solar calculations.
      ! -------------------------------------------------------
      iN = 0
      INPUT_VARS_2: do k=1,NumInp
         if (SlicesInp(k) == 0) cycle

         i1 = iN + 1

         if (NamesInp(k)=="AERO") then

            _ASSERT(size(AEROSOL_EXT,3)==LM,'mal-dimensioned AEROSOL_EXT')
            _ASSERT(size(AEROSOL_SSA,3)==LM,'mal-dimensioned AEROSOL_SSA')
            _ASSERT(size(AEROSOL_ASY,3)==LM,'mal-dimensioned AEROSOL_ASY')

            allocate(BUF_AEROSOL(size(AEROSOL_EXT,1), &
                                 size(AEROSOL_EXT,2), &
                                 size(AEROSOL_EXT,3)), __STAT__)

            ! pack extinctions
            BUF_AEROSOL = MAPL_UNDEF
            do j=1,NUM_BANDS_SOLAR
               BUF_AEROSOL = AEROSOL_EXT(:,:,:,j)
               call PackIt(BufInp(i1+(j-1)*LM*NumMax),BUF_AEROSOL,daytime,NumMax,HorzDims,LM)
            end do
            iN = i1 + NumMax*LM*NUM_BANDS_SOLAR - 1
            ptr3(1:NumMax,1:LM,1:NUM_BANDS_SOLAR) => BufInp(i1:iN)
            BUFIMP_AEROSOL_EXT => ptr3(1:Num2do,:,:)

            ! pack single scattering albedos
            i1 = iN + 1
            BUF_AEROSOL = MAPL_UNDEF
            do j=1,NUM_BANDS_SOLAR
               BUF_AEROSOL = AEROSOL_SSA(:,:,:,j)
               call PackIt(BufInp(i1+(j-1)*LM*NumMax),BUF_AEROSOL,daytime,NumMax,HorzDims,LM)
            end do
            iN = i1 + NumMax*LM*NUM_BANDS_SOLAR - 1
            ptr3(1:NumMax,1:LM,1:NUM_BANDS_SOLAR) => BufInp(i1:iN)
            BUFIMP_AEROSOL_SSA => ptr3(1:Num2do,:,:)

            ! pack asymmetry factors
            i1 = iN + 1
            BUF_AEROSOL = MAPL_UNDEF
            do j=1,NUM_BANDS_SOLAR
               BUF_AEROSOL = AEROSOL_ASY(:,:,:,j)
               call PackIt(BufInp(i1+(j-1)*LM*NumMax),BUF_AEROSOL,daytime,NumMax,HorzDims,LM)
            end do
            iN = i1 + NumMax*LM*NUM_BANDS_SOLAR - 1
            ptr3(1:NumMax,1:LM,1:NUM_BANDS_SOLAR) => BufInp(i1:iN)
            BUFIMP_AEROSOL_ASY => ptr3(1:Num2do,:,:)

            deallocate(BUF_AEROSOL, __STAT__)

         else  ! Non-aerosol imports

            if (SlicesInp(k) /= 1) then

               ! pack 3D imports
               call ESMFL_StateGetPointerToData(IMPORT,ptr3,NamesInp(k),__RC__)
               call PackIt(BufInp(i1),ptr3,daytime,NumMax,HorzDims,size(ptr3,3))
               iN = i1 + NumMax*size(ptr3,3) - 1

            else  ! case(MAPL_DIMSHORZONLY)

               ! pack auxilliary variables
               if (NamesInp(k) == 'Ig') then
                  call PackIt(BufInp(i1),real(Ig),daytime,NumMax,HorzDims,1)
               else if (NamesInp(k) == 'Jg') then
                  call PackIt(BufInp(i1),real(Jg),daytime,NumMax,HorzDims,1)
               else if (NamesInp(k) == 'LATS') then
                  call PackIt(BufInp(i1),LATS,    daytime,NumMax,HorzDims,1)
               else if (NamesInp(k) == 'SLR') then
                  call PackIt(BufInp(i1),SLR,     daytime,NumMax,HorzDims,1)
               else if (NamesInp(k) == 'ZTH') then
                  call PackIt(BufInp(i1),ZTH,     daytime,NumMax,HorzDims,1)
               else
                  ! pack 2D imports
                  call ESMFL_StateGetPointerToData(IMPORT,ptr2,NamesInp(k),__RC__)
                  call PackIt(BufInp(i1),ptr2,daytime,NumMax,HorzDims,1)
               end if
               iN = i1 + NumMax - 1

            end if

            ! Handles for the working input (Import) variables.
            ! These use Fortran 2003 syntax for reshaping a 1D
            ! vector into a higher rank array.
            !--------------------------------------------------
            ptr2(1:NumMax,1:SlicesInp(k)) => BufInp(i1:iN)

            select case(NamesInp(k))
               case('PLE')
                  PLE   => ptr2(1:Num2do,:)
               case('TS')
                  TS    => ptr2(1:Num2do,1)
               case('CH4')
                  CH4   => ptr2(1:Num2do,:)
               case('N2O')
                  N2O   => ptr2(1:Num2do,:)
               case('T')
                  T     => ptr2(1:Num2do,:)
               case('QV')
                  Q     => ptr2(1:Num2do,:)
               case('OX')
                  OX    => ptr2(1:Num2do,:)
               case('FCLD')
                  CL    => ptr2(1:Num2do,:)
               case('QL')
                  QL    => ptr2(1:Num2do,:)
               case('QI')
                  QI    => ptr2(1:Num2do,:)
               case('QR')
                  QR    => ptr2(1:Num2do,:)
               case('QS')
                  QS    => ptr2(1:Num2do,:)
               case('RL')
                  RL    => ptr2(1:Num2do,:)
               case('RI')
                  RI    => ptr2(1:Num2do,:)
               case('RR')
                  RR    => ptr2(1:Num2do,:)
               case('RS')
                  RS    => ptr2(1:Num2do,:)
               case('ALBVR')
                  ALBVR => ptr2(1:Num2do,1)
               case('ALBVF')
                  ALBVF => ptr2(1:Num2do,1)
               case('ALBNR')
                  ALBNR => ptr2(1:Num2do,1)
               case('ALBNF')
                  ALBNF => ptr2(1:Num2do,1)
               case('Ig')
                  Ig1D  => ptr2(1:Num2do,1)
               case('Jg')
                  Jg1D  => ptr2(1:Num2do,1)
               case('LATS')
                  ALAT  => ptr2(1:Num2do,1)
               case('SLR')
                  SLR1D => ptr2(1:Num2do,1)
               case('ZTH')
                  ZT    => ptr2(1:Num2do,1)
            end select

         end if

      enddo INPUT_VARS_2

      ! Load balance the Inputs
      ! -----------------------

      call MAPL_TimerOn(MAPL,"--DISTRIBUTE")
      if (LoadBalance) call MAPL_BalanceWork(BufInp,NumMax,Direction=MAPL_Distribute,Handle=SolarBalanceHandle,__RC__)
      call MAPL_TimerOff(MAPL,"--DISTRIBUTE")

! @@@@@@@@@@@@@@@@@@@@@@
! @@@ InOuts/Outputs @@@
! @@@@@@@@@@@@@@@@@@@@@@

      ! Count the 2D slices in the Int (InOut/Out) vars and calc
      ! the required length of their 1D buffers (BufInOut/BufOut).
      ! ----------------------------------------------------------
      INT_VARS_1: do k=1,NumInt

         ! InOut or Out?
         call MAPL_VarSpecGet(InternalSpec(k), &
            SHORT_NAME=short_name, DIMS=dims, UNGRIDDED_DIMS=ugdims, __RC__)
         ! later FAR variables will be InOut ... for now there are no InOut vars
         IntInOut(k) = .false.

         ! save properties
         NamesInt(k) = short_name
         rgDim(k) = dims

         ! Skip vertical only variables. They dont require
         ! load-balancing since they have no horizontal dimension.
         if (dims == MAPL_DIMSVERTONLY) then
            SlicesInt(k) = 0
            cycle
         end if

         ! Exclude unused internals
         if (.not. include_aerosols .and.                                    &
               ('FSWN'       == short_name .or.    'FSCN' == short_name .or. &
                'FSWUN'      == short_name .or.   'FSCUN' == short_name .or. &
                'FSWBANDN'   == short_name)                                  &
            .or. include_aerosols .and.                                      &
               ('FSWNAN'     == short_name .or.  'FSCNAN' == short_name .or. &
                'FSWUNAN'    == short_name .or. 'FSCUNAN' == short_name .or. &
                'FSWBANDNAN' == short_name)                                  &
            ) then
            SlicesInt(k) = 0
            cycle
         end if

         ! Don't calculate D[RF]BANDN for .not. include_aerosols
         if (.not. include_aerosols) then
            if (short_name == 'DRBANDN' .or. short_name == 'DFBANDN') then
               SlicesInt(k) = 0
               cycle
            end if
         end if

         if (associated(ugdims)) then
            ! ungridded dims are present, make sure just one
            _ASSERT(size(ugdims)==1,'Only one ungridded dimension allowed')
            ugDim(k) = ugdims(1)
            select case(dims)
               case(MAPL_DIMSHORZVERT)
                  call ESMFL_StateGetPointerToData(INTERNAL,ptr4,NamesInt(k),__RC__)
                  SlicesInt(k) = size(ptr4,3) * ugDim(k)
               case(MAPL_DIMSHORZONLY)
                  SlicesInt(k) = ugDim(k)
               case default
                  _FAIL('invalid dimension for Internal')
            end select
         else
            ! no ungridded dimension
            ugDim(k) = 0
            select case(dims)
               case(MAPL_DIMSHORZVERT)
                  call ESMFL_StateGetPointerToData(INTERNAL,ptr3,NamesInt(k),__RC__)
                  SlicesInt(k) = size(ptr3,3)
               case(MAPL_DIMSHORZONLY)
                  SlicesInt(k) = 1
               case default
                  _FAIL('invalid dimension for Internal')
            end select
         end if

      enddo INT_VARS_1

      ! Allocate buffers with enough space to hold both the unbalanced
      ! and balanced data on the local PE for InOut/Out vars
      ! --------------------------------------------------------------
      allocate(BufInOut(NumMax*sum(SlicesInt,MASK=IntInOut)),__STAT__)
      BufInOut = MAPL_UNDEF
      allocate(BufOut(NumMax*sum(SlicesInt,MASK=.not.IntInOut)),__STAT__)
      BufOut = MAPL_UNDEF

      ! Loop over Internals (InOuts/Outs), packing them into buffers
      ! that will be load balanced and used in the solar calculations.
      ! --------------------------------------------------------------
      iNInOut = 0; iNOut = 0
      INT_VARS_2: do k=1,NumInt
         if (SlicesInt(k) == 0) cycle

         if (IntInOut(k)) then
            buf => bufInOut; pi1 => i1InOut; piN => iNInOut
         else
            buf => bufOut;   pi1 => i1Out;   piN => iNOut
         endif
         pi1 = piN + 1

         if (ugDim(k) > 0) then  ! has ungridded dimensions

            select case(rgDim(k))
               case(MAPL_DIMSHORZVERT)
                  call ESMFL_StateGetPointerToData(INTERNAL,ptr4,NamesInt(k),__RC__)
                  do j=1,ugDim(k)
!pmn compiler       call PackIt(Buf(pi1+(j-1)*size(ptr4,3)*NumMax),ptr4(:,:,:,j),daytime,NumMax,HorzDims,size(ptr4,3))
                    if (IntInOut(k)) then
                      call PackIt(BufInOut(pi1+(j-1)*size(ptr4,3)*NumMax),ptr4(:,:,:,j),daytime,NumMax,HorzDims,size(ptr4,3))
                    else
                      call PackIt(BufOut  (pi1+(j-1)*size(ptr4,3)*NumMax),ptr4(:,:,:,j),daytime,NumMax,HorzDims,size(ptr4,3))
                    endif
                  end do
                  piN = pi1 + NumMax*size(ptr4,3)*ugDim(k) - 1
                  ptr3(1:NumMax,1:size(ptr4,3),1:ugDim(k)) => Buf(pi1:piN)
               case(MAPL_DIMSHORZONLY)
                  call ESMFL_StateGetPointerToData(INTERNAL,ptr3,NamesInt(k),__RC__)
!pmn compiler     call PackIt(Buf(pi1),ptr3,daytime,NumMax,HorzDims,ugDim(k))
                  if (IntInOut(k)) then
                    call PackIt(BufInOut(pi1),ptr3,daytime,NumMax,HorzDims,ugDim(k))
                  else
                    call PackIt(BufOut  (pi1),ptr3,daytime,NumMax,HorzDims,ugDim(k))
                  endif
                  piN = pi1 + NumMax*ugDim(k) - 1
                  ptr2(1:NumMax,1:ugDim(k)) => Buf(pi1:piN)
            end select

         else  ! no ungridded dimensions

            select case(rgDim(k))
               case(MAPL_DIMSHORZVERT)
                  call ESMFL_StateGetPointerToData(INTERNAL,ptr3,NamesInt(k),__RC__)
!pmn compiler     call PackIt(Buf(pi1),ptr3,daytime,NumMax,HorzDims,size(ptr3,3))
                  if (IntInOut(k)) then
                     call PackIt(BufInOut(pi1),ptr3,daytime,NumMax,HorzDims,size(ptr3,3))
                  else
                     call PackIt(BufOut  (pi1),ptr3,daytime,NumMax,HorzDims,size(ptr3,3))
                  endif
                  piN = pi1 + NumMax*size(ptr3,3) - 1
               case(MAPL_DIMSHORZONLY)
                  call ESMFL_StateGetPointerToData(INTERNAL,ptr2,NamesInt(k),__RC__)
!pmn compiler     call PackIt(Buf(pi1),ptr2,daytime,NumMax,HorzDims,1)
                  if (IntInOut(k)) then
                     call PackIt(BufInOut(pi1),ptr2,daytime,NumMax,HorzDims,1)
                  else
                     call PackIt(BufOut  (pi1),ptr2,daytime,NumMax,HorzDims,1)
                  endif
                  piN = pi1 + NumMax - 1
            end select
            ptr2(1:NumMax,1:SlicesInt(k)) => Buf(pi1:piN)

         end if

         ! Handles for the working InOut/Out variables.
         ! These have an inner dimension of the balanced work.
         ! ---------------------------------------------------
         select case(NamesInt(k))
            case('FSWN')
               FSW       => ptr2(1:Num2do,:)
            case('FSCN')
               FSC       => ptr2(1:Num2do,:)
            case('FSWUN')
               FSWU      => ptr2(1:Num2do,:)
            case('FSCUN')
               FSCU      => ptr2(1:Num2do,:)
            case('FSWBANDN')
               FSWBAND   => ptr2(1:Num2do,:)
            case('DRUVRN')
               UVRR      => ptr2(1:Num2do,1)
            case('DFUVRN')
               UVRF      => ptr2(1:Num2do,1)
            case('DRPARN')
               PARR      => ptr2(1:Num2do,1)
            case('DFPARN')
               PARF      => ptr2(1:Num2do,1)
            case('DRNIRN')
               NIRR      => ptr2(1:Num2do,1)
            case('DFNIRN')
               NIRF      => ptr2(1:Num2do,1)
            case('DRBANDN')
               DRBAND    => ptr2(1:Num2do,:)
            case('DFBANDN')
               DFBAND    => ptr2(1:Num2do,:)
            case('FSWNAN')
               FSWA      => ptr2(1:Num2do,:)
            case('FSCNAN')
               FSCA      => ptr2(1:Num2do,:)
            case('FSWUNAN')
               FSWUA     => ptr2(1:Num2do,:)
            case('FSCUNAN')
               FSCUA     => ptr2(1:Num2do,:)
            case('FSWBANDNAN')
               FSWBANDA  => ptr2(1:Num2do,:)
            case('COSZSW')
               COSZSW    => ptr2(1:Num2do,1)
            case('CLDTTSW')
               CLDTS     => ptr2(1:Num2do,1)
            case('CLDHISW')
               CLDHS     => ptr2(1:Num2do,1)
            case('CLDMDSW')
               CLDMS     => ptr2(1:Num2do,1)
            case('CLDLOSW')
               CLDLS     => ptr2(1:Num2do,1)
#ifdef SOLAR_RADVAL
            case('TAUTTPAR')
               TAUTP     => ptr2(1:Num2do,1)
            case('TAUHIPAR')
               TAUHP     => ptr2(1:Num2do,1)
            case('TAUMDPAR')
               TAUMP     => ptr2(1:Num2do,1)
            case('TAULOPAR')
               TAULP     => ptr2(1:Num2do,1)
#endif
            case('COTTTPAR')
               COTTP     => ptr2(1:Num2do,1)
            case('COTHIPAR')
               COTHP     => ptr2(1:Num2do,1)
            case('COTMDPAR')
               COTMP     => ptr2(1:Num2do,1)
            case('COTLOPAR')
               COTLP     => ptr2(1:Num2do,1)
            case('COTDENTTPAR')
               COTDTP    => ptr2(1:Num2do,1)
            case('COTDENHIPAR')
               COTDHP    => ptr2(1:Num2do,1)
            case('COTDENMDPAR')
               COTDMP    => ptr2(1:Num2do,1)
            case('COTDENLOPAR')
               COTDLP    => ptr2(1:Num2do,1)
            case('COTNUMTTPAR')
               COTNTP    => ptr2(1:Num2do,1)
            case('COTNUMHIPAR')
               COTNHP    => ptr2(1:Num2do,1)
            case('COTNUMMDPAR')
               COTNMP    => ptr2(1:Num2do,1)
            case('COTNUMLOPAR')
               COTNLP    => ptr2(1:Num2do,1)
#ifdef SOLAR_RADVAL
            case('COTDSDENTTPAR')
               CDSDTP    => ptr2(1:Num2do,1)
            case('COTDSDENHIPAR')
               CDSDHP    => ptr2(1:Num2do,1)
            case('COTDSDENMDPAR')
               CDSDMP    => ptr2(1:Num2do,1)
            case('COTDSDENLOPAR')
               CDSDLP    => ptr2(1:Num2do,1)
            case('COTDSNUMTTPAR')
               CDSNTP    => ptr2(1:Num2do,1)
            case('COTDSNUMHIPAR')
               CDSNHP    => ptr2(1:Num2do,1)
            case('COTDSNUMMDPAR')
               CDSNMP    => ptr2(1:Num2do,1)
            case('COTDSNUMLOPAR')
               CDSNLP    => ptr2(1:Num2do,1)
            case('COTLDENTTPAR')
               COTLDTP    => ptr2(1:Num2do,1)
            case('COTLDENHIPAR')
               COTLDHP    => ptr2(1:Num2do,1)
            case('COTLDENMDPAR')
               COTLDMP    => ptr2(1:Num2do,1)
            case('COTLDENLOPAR')
               COTLDLP    => ptr2(1:Num2do,1)
            case('COTLNUMTTPAR')
               COTLNTP    => ptr2(1:Num2do,1)
            case('COTLNUMHIPAR')
               COTLNHP    => ptr2(1:Num2do,1)
            case('COTLNUMMDPAR')
               COTLNMP    => ptr2(1:Num2do,1)
            case('COTLNUMLOPAR')
               COTLNLP    => ptr2(1:Num2do,1)
            case('COTLDSDENTTPAR')
               CDSLDTP    => ptr2(1:Num2do,1)
            case('COTLDSDENHIPAR')
               CDSLDHP    => ptr2(1:Num2do,1)
            case('COTLDSDENMDPAR')
               CDSLDMP    => ptr2(1:Num2do,1)
            case('COTLDSDENLOPAR')
               CDSLDLP    => ptr2(1:Num2do,1)
            case('COTLDSNUMTTPAR')
               CDSLNTP    => ptr2(1:Num2do,1)
            case('COTLDSNUMHIPAR')
               CDSLNHP    => ptr2(1:Num2do,1)
            case('COTLDSNUMMDPAR')
               CDSLNMP    => ptr2(1:Num2do,1)
            case('COTLDSNUMLOPAR')
               CDSLNLP    => ptr2(1:Num2do,1)
            case('COTIDENTTPAR')
               COTIDTP    => ptr2(1:Num2do,1)
            case('COTIDENHIPAR')
               COTIDHP    => ptr2(1:Num2do,1)
            case('COTIDENMDPAR')
               COTIDMP    => ptr2(1:Num2do,1)
            case('COTIDENLOPAR')
               COTIDLP    => ptr2(1:Num2do,1)
            case('COTINUMTTPAR')
               COTINTP    => ptr2(1:Num2do,1)
            case('COTINUMHIPAR')
               COTINHP    => ptr2(1:Num2do,1)
            case('COTINUMMDPAR')
               COTINMP    => ptr2(1:Num2do,1)
            case('COTINUMLOPAR')
               COTINLP    => ptr2(1:Num2do,1)
            case('COTIDSDENTTPAR')
               CDSIDTP    => ptr2(1:Num2do,1)
            case('COTIDSDENHIPAR')
               CDSIDHP    => ptr2(1:Num2do,1)
            case('COTIDSDENMDPAR')
               CDSIDMP    => ptr2(1:Num2do,1)
            case('COTIDSDENLOPAR')
               CDSIDLP    => ptr2(1:Num2do,1)
            case('COTIDSNUMTTPAR')
               CDSINTP    => ptr2(1:Num2do,1)
            case('COTIDSNUMHIPAR')
               CDSINHP    => ptr2(1:Num2do,1)
            case('COTIDSNUMMDPAR')
               CDSINMP    => ptr2(1:Num2do,1)
            case('COTIDSNUMLOPAR')
               CDSINLP    => ptr2(1:Num2do,1)
            case('SSALDENLOPAR')
               SSALDLP    => ptr2(1:Num2do,1)
            case('SSALNUMLOPAR')
               SSALNLP    => ptr2(1:Num2do,1)
            case('SSAIDENLOPAR')
               SSAIDLP    => ptr2(1:Num2do,1)
            case('SSAINUMLOPAR')
               SSAINLP    => ptr2(1:Num2do,1)
            case('ASMLDENLOPAR')
               ASMLDLP    => ptr2(1:Num2do,1)
            case('ASMLNUMLOPAR')
               ASMLNLP    => ptr2(1:Num2do,1)
            case('ASMIDENLOPAR')
               ASMIDLP    => ptr2(1:Num2do,1)
            case('ASMINUMLOPAR')
               ASMINLP    => ptr2(1:Num2do,1)
            case('SSALDENMDPAR')
               SSALDMP    => ptr2(1:Num2do,1)
            case('SSALNUMMDPAR')
               SSALNMP    => ptr2(1:Num2do,1)
            case('SSAIDENMDPAR')
               SSAIDMP    => ptr2(1:Num2do,1)
            case('SSAINUMMDPAR')
               SSAINMP    => ptr2(1:Num2do,1)
            case('ASMLDENMDPAR')
               ASMLDMP    => ptr2(1:Num2do,1)
            case('ASMLNUMMDPAR')
               ASMLNMP    => ptr2(1:Num2do,1)
            case('ASMIDENMDPAR')
               ASMIDMP    => ptr2(1:Num2do,1)
            case('ASMINUMMDPAR')
               ASMINMP    => ptr2(1:Num2do,1)
            case('SSALDENHIPAR')
               SSALDHP    => ptr2(1:Num2do,1)
            case('SSALNUMHIPAR')
               SSALNHP    => ptr2(1:Num2do,1)
            case('SSAIDENHIPAR')
               SSAIDHP    => ptr2(1:Num2do,1)
            case('SSAINUMHIPAR')
               SSAINHP    => ptr2(1:Num2do,1)
            case('ASMLDENHIPAR')
               ASMLDHP    => ptr2(1:Num2do,1)
            case('ASMLNUMHIPAR')
               ASMLNHP    => ptr2(1:Num2do,1)
            case('ASMIDENHIPAR')
               ASMIDHP    => ptr2(1:Num2do,1)
            case('ASMINUMHIPAR')
               ASMINHP    => ptr2(1:Num2do,1)
            case('SSALDENTTPAR')
               SSALDTP    => ptr2(1:Num2do,1)
            case('SSALNUMTTPAR')
               SSALNTP    => ptr2(1:Num2do,1)
            case('SSAIDENTTPAR')
               SSAIDTP    => ptr2(1:Num2do,1)
            case('SSAINUMTTPAR')
               SSAINTP    => ptr2(1:Num2do,1)
            case('ASMLDENTTPAR')
               ASMLDTP    => ptr2(1:Num2do,1)
            case('ASMLNUMTTPAR')
               ASMLNTP    => ptr2(1:Num2do,1)
            case('ASMIDENTTPAR')
               ASMIDTP    => ptr2(1:Num2do,1)
            case('ASMINUMTTPAR')
               ASMINTP    => ptr2(1:Num2do,1)
            case('SSALDSDENLOPAR')
               SDSLDLP    => ptr2(1:Num2do,1)
            case('SSALDSNUMLOPAR')
               SDSLNLP    => ptr2(1:Num2do,1)
            case('SSAIDSDENLOPAR')
               SDSIDLP    => ptr2(1:Num2do,1)
            case('SSAIDSNUMLOPAR')
               SDSINLP    => ptr2(1:Num2do,1)
            case('ASMLDSDENLOPAR')
               ADSLDLP    => ptr2(1:Num2do,1)
            case('ASMLDSNUMLOPAR')
               ADSLNLP    => ptr2(1:Num2do,1)
            case('ASMIDSDENLOPAR')
               ADSIDLP    => ptr2(1:Num2do,1)
            case('ASMIDSNUMLOPAR')
               ADSINLP    => ptr2(1:Num2do,1)
            case('SSALDSDENMDPAR')
               SDSLDMP    => ptr2(1:Num2do,1)
            case('SSALDSNUMMDPAR')
               SDSLNMP    => ptr2(1:Num2do,1)
            case('SSAIDSDENMDPAR')
               SDSIDMP    => ptr2(1:Num2do,1)
            case('SSAIDSNUMMDPAR')
               SDSINMP    => ptr2(1:Num2do,1)
            case('ASMLDSDENMDPAR')
               ADSLDMP    => ptr2(1:Num2do,1)
            case('ASMLDSNUMMDPAR')
               ADSLNMP    => ptr2(1:Num2do,1)
            case('ASMIDSDENMDPAR')
               ADSIDMP    => ptr2(1:Num2do,1)
            case('ASMIDSNUMMDPAR')
               ADSINMP    => ptr2(1:Num2do,1)
            case('SSALDSDENHIPAR')
               SDSLDHP    => ptr2(1:Num2do,1)
            case('SSALDSNUMHIPAR')
               SDSLNHP    => ptr2(1:Num2do,1)
            case('SSAIDSDENHIPAR')
               SDSIDHP    => ptr2(1:Num2do,1)
            case('SSAIDSNUMHIPAR')
               SDSINHP    => ptr2(1:Num2do,1)
            case('ASMLDSDENHIPAR')
               ADSLDHP    => ptr2(1:Num2do,1)
            case('ASMLDSNUMHIPAR')
               ADSLNHP    => ptr2(1:Num2do,1)
            case('ASMIDSDENHIPAR')
               ADSIDHP    => ptr2(1:Num2do,1)
            case('ASMIDSNUMHIPAR')
               ADSINHP    => ptr2(1:Num2do,1)
            case('SSALDSDENTTPAR')
               SDSLDTP    => ptr2(1:Num2do,1)
            case('SSALDSNUMTTPAR')
               SDSLNTP    => ptr2(1:Num2do,1)
            case('SSAIDSDENTTPAR')
               SDSIDTP    => ptr2(1:Num2do,1)
            case('SSAIDSNUMTTPAR')
               SDSINTP    => ptr2(1:Num2do,1)
            case('ASMLDSDENTTPAR')
               ADSLDTP    => ptr2(1:Num2do,1)
            case('ASMLDSNUMTTPAR')
               ADSLNTP    => ptr2(1:Num2do,1)
            case('ASMIDSDENTTPAR')
               ADSIDTP    => ptr2(1:Num2do,1)
            case('ASMIDSNUMTTPAR')
               ADSINTP    => ptr2(1:Num2do,1)
            case('FORLDENLOPAR')
               FORLDLP    => ptr2(1:Num2do,1)
            case('FORLNUMLOPAR')
               FORLNLP    => ptr2(1:Num2do,1)
            case('FORIDENLOPAR')
               FORIDLP    => ptr2(1:Num2do,1)
            case('FORINUMLOPAR')
               FORINLP    => ptr2(1:Num2do,1)
            case('FORLDENMDPAR')
               FORLDMP    => ptr2(1:Num2do,1)
            case('FORLNUMMDPAR')
               FORLNMP    => ptr2(1:Num2do,1)
            case('FORIDENMDPAR')
               FORIDMP    => ptr2(1:Num2do,1)
            case('FORINUMMDPAR')
               FORINMP    => ptr2(1:Num2do,1)
            case('FORLDENHIPAR')
               FORLDHP    => ptr2(1:Num2do,1)
            case('FORLNUMHIPAR')
               FORLNHP    => ptr2(1:Num2do,1)
            case('FORIDENHIPAR')
               FORIDHP    => ptr2(1:Num2do,1)
            case('FORINUMHIPAR')
               FORINHP    => ptr2(1:Num2do,1)
            case('FORLDENTTPAR')
               FORLDTP    => ptr2(1:Num2do,1)
            case('FORLNUMTTPAR')
               FORLNTP    => ptr2(1:Num2do,1)
            case('FORIDENTTPAR')
               FORIDTP    => ptr2(1:Num2do,1)
            case('FORINUMTTPAR')
               FORINTP    => ptr2(1:Num2do,1)
#endif
         end select

      enddo INT_VARS_2

      ! Load balance the InOuts for Input
      !----------------------------------
      call MAPL_TimerOn(MAPL,"--DISTRIBUTE")
      if (size(BufInOut) > 0) then
         if (LoadBalance) call MAPL_BalanceWork(BufInOut,NumMax,Direction=MAPL_Distribute,Handle=SolarBalanceHandle,__RC__)
      end if
      call MAPL_TimerOff(MAPL,"--DISTRIBUTE")

      call MAPL_TimerOff(MAPL,"-BALANCE")

! Do shortwave calculations on a list of soundings
!-------------------------------------------------

      call MAPL_TimerOn(MAPL,"-MISC")

      ! report cosine solar zenith angle actually used by REFRESH
      COSZSW = ZT

      ! save soon-to-be-calculated fluxes to correct set of internals
      if (.not. include_aerosols) then
         FSW     => FSWA
         FSC     => FSCA
         FSWU    => FSWUA
         FSCU    => FSCUA
         FSWBAND => FSWBANDA
      end if

      ! Option to force binary clouds for SW
      call MAPL_GetResource(MAPL,ibinary,"RADSW_BINARY_CLOUDS:",DEFAULT=0,__RC__)
      if (ibinary /= 0) where (CL > 0.) CL = 1.

      ! Prepare auxilliary variables
      ! ----------------------------

      allocate(RH(size(Q,1),size(Q,2)),__STAT__)
      allocate(PL(size(Q,1),size(Q,2)),__STAT__)
      allocate(PLhPa(size(PLE,1),size(PLE,2)),__STAT__)

      PL = 0.5*(PLE(:,:UBOUND(PLE,2)-1)+PLE(:,LBOUND(PLE,2)+1:))
      RH = Q/MAPL_EQSAT(T,PL=PL)
      PLhPa = PLE * 0.01

      ! Water amounts and effective radii are in arrays indexed by species
      !-------------------------------------------------------------------

      allocate(QQ3 (size(Q,1),size(Q,2),4),__STAT__)
      allocate(RR3 (size(Q,1),size(Q,2),4),__STAT__)

      ! In-cloud water contents
      QQ3(:,:,1) = QI
      QQ3(:,:,2) = QL
      QQ3(:,:,3) = QR
      QQ3(:,:,4) = QS

      ! Effective radii [microns]
      WHERE (RI == MAPL_UNDEF) RI = 36.e-6
      WHERE (RL == MAPL_UNDEF) RL = 14.e-6
      WHERE (RR == MAPL_UNDEF) RR = 50.e-6
      WHERE (RS == MAPL_UNDEF) RS = 50.e-6
      RR3(:,:,1) = RI*1.e6
      RR3(:,:,2) = RL*1.e6
      RR3(:,:,3) = RR*1.e6
      RR3(:,:,4) = RS*1.e6

      ! Convert odd oxygen, which is the model prognostic, to ozone
      !------------------------------------------------------------

      allocate(O3 (size(Q,1),size(Q,2)),__STAT__)

      O3 = OX
      WHERE(PL < 100.)
         O3 = O3 * EXP(-1.5*(LOG10(PL)-2.)**2)
      ENDWHERE

      ! SORAD expects non-negative ozone fraction by MASS
      !--------------------------------------------------

      O3 = O3 * (MAPL_O3MW / MAPL_AIRMW)
      O3 = MAX(O3, 0.00)

      ! ------------------
      ! Begin aerosol code
      ! ------------------

      allocate(TAUA(size(Q,1),size(Q,2),NUM_BANDS_SOLAR),__STAT__)
      allocate(SSAA(size(Q,1),size(Q,2),NUM_BANDS_SOLAR),__STAT__)
      allocate(ASYA(size(Q,1),size(Q,2),NUM_BANDS_SOLAR),__STAT__)

      ! Zero out aerosol arrays.
      ! If num_aero_vars == 0, these zeroes are used inside code.
      TAUA = 0.
      SSAA = 0.
      ASYA = 0.

      ! If we have aerosols, load them.
      if (num_aero_vars > 0) then
         TAUA = BUFIMP_AEROSOL_EXT
         SSAA = BUFIMP_AEROSOL_SSA
         ASYA = BUFIMP_AEROSOL_ASY
      end if

      call MAPL_TimerOff(MAPL,"-MISC")

      ! Call the requested Shortwave scheme
      ! -----------------------------------

   SCHEME: if (USE_CHOU) then
      call shrtwave(                                          &
                     PLhPa, T, Q, O3, CO2, ZT,                &
                     QQ3, RR3, CL,                            &
                     LCLDMH,LCLDLM,                           &
                     ALBVR, ALBVF, ALBNR, ALBNF,              &
                     TAUA,SSAA,ASYA,                          &

                     FSW ,FSC ,                               &
                     NIRR,NIRF,PARR,PARF,UVRR,UVRF,           &
                     FSWU,FSCU,                               &
                     FSWBAND,                                 &
                     SOLAR_TO_OBIO .and. include_aerosols,    &
                     DRBAND, DFBAND,                          &
                     __RC__                                   )

   else if (USE_RRTMGP) then

! helper for testing RRTMGP error status on return
! allows line number reporting cf. original call method
#define TEST_(A) error_msg = A; if (trim(error_msg)/="") then; _FAIL("RRTMGP Error: "//trim(error_msg)); endif

      call MAPL_TimerOn(MAPL,"-RRTMGP",__RC__)

      ! number of columns after load balancing
      ncol = size(Q,1)

      ! absorbing gas names
      error_msg = gas_concs%init([character(3) :: &
        'h2o','co2','o3','n2o','co','ch4','o2','n2'])
      TEST_(error_msg)

      ! load gas concentrations (volume mixing ratios)
      ! "constant" gases
      TEST_(gas_concs%set_vmr('n2' , real(N2 ,kind=wp)))
      TEST_(gas_concs%set_vmr('o2' , real(O2 ,kind=wp)))
      TEST_(gas_concs%set_vmr('co2', real(CO2,kind=wp)))
      TEST_(gas_concs%set_vmr('co' , real(CO ,kind=wp)))
      ! variable gases
      ! (ozone converted from mass mixing ratio, water vapor from specific humidity)
      TEST_(gas_concs%set_vmr('ch4', real(CH4                             ,kind=wp)))
      TEST_(gas_concs%set_vmr('n2o', real(N2O                             ,kind=wp)))
      TEST_(gas_concs%set_vmr('o3' , real(O3      *(MAPL_AIRMW/MAPL_O3MW ),kind=wp)))
      TEST_(gas_concs%set_vmr('h2o', real(Q/(1.-Q)*(MAPL_AIRMW/MAPL_H2OMW),kind=wp)))

      ! access RRTMGP internal state from the GC
      call ESMF_UserCompGetInternalState(GC, 'RRTMGP_state', wrap, status)
      VERIFY_(status)
      rrtmgp_state => wrap%ptr

      ! initialize k-distribution if not already done
      call MAPL_GetResource( &
        MAPL, k_dist_file, "RRTMGP_GAS_SW:", &
        DEFAULT='rrtmgp-gas-sw-g112.nc',__RC__)
      if (.not. rrtmgp_state%initialized) then
        ! gas_concs needed only to access required gas names
        call MAPL_TimerOn(MAPL,"--RRTMGP_IO_GAS",__RC__)
        call load_and_init(rrtmgp_state%k_dist,trim(k_dist_file),gas_concs)
        call MAPL_TimerOff(MAPL,"--RRTMGP_IO_GAS",__RC__)
        if (.not. rrtmgp_state%k_dist%source_is_external()) then
          TEST_('RRTMGP-SW: does not seem to be SW')
        endif
        rrtmgp_state%initialized = .true.
      endif

      ! access by shorter name
      k_dist => rrtmgp_state%k_dist

      ! adjust sun for current faculae and sunspots (tsi scaling done later)
      ! for the moment we are requiring MG and SB indicies from the NRLSSI2 file
      _ASSERT(SolCycFileName /= '/dev/null' .and. USE_NRLSSI2, 'RRTMGP-SW: MG and SB not available')
      error_msg = k_dist%set_solar_variability(real(MG,kind=wp),real(SB,kind=wp))
      TEST_(error_msg)

      ! spectral dimensions
      ngpt = k_dist%get_ngpt()
      nbnd = k_dist%get_nband()
      _ASSERT(nbnd == NB_RRTMGP, 'RRTMGP-SW: expected different number of bands')

      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! For reference, comparison of RRTMG and RRTMGP bands:
      ! from RRTMG:
      ! wavenum1(:) = (/2600., 3250., 4000., 4650., 5150., 6150., 7700.,  8050.,12850.,16000.,22650.,29000.,38000.,  820./)
      ! wavenum2(:) = (/3250., 4000., 4650., 5150., 6150., 7700., 8050., 12850.,16000.,22650.,29000.,38000.,50000., 2600./)
      ! from RRTMGP:
      ! write(*,*) 'band_lims_wvn(2,nbnd):', k_dist%get_band_lims_wavenumber() ! reordered below
      !               820., 2680., 3250., 4000., 4650., 5150., 6150., 7700.,  8050., 12850., 16000., 22650., 29000., 38000.
      !              2680., 3250., 4000., 4650., 5150., 6150., 7700., 8050., 12850., 16000., 22650., 29000., 38000., 50000.
      ! clearly there are some differences ... so aerosol tables were redone
      !   mainly band 14 becomes band 1, plus small change in wavenumber upper limit of that band only
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! gpoint limits for each band
      allocate (band_lims_gpt(2,nbnd),__STAT__)
      band_lims_gpt = k_dist%get_band_lims_gpoint()

      ! dummy array (see later)
      allocate(dummy_wp(ncol,LM),source=0._wp, __STAT__)

      ! allocate input arrays
      allocate(tsi(ncol), mu0(ncol), __STAT__)
      allocate(sfc_alb_dir(nbnd,ncol), sfc_alb_dif(nbnd,ncol), __STAT__)
      allocate(p_lay(ncol,LM), t_lay(ncol,LM), dp_wp(ncol,LM), __STAT__)
      allocate(dzmid(ncol,LM-1), __STAT__)
      allocate(p_lev(ncol,LM+1), __STAT__)

      ! load input arrays ...

      ! solar inputs:
      ! 1. cosine of solar zenith angle, mu0
      ! Note: this ZT, ultimately from ZTH of MAPL_SunGetInsolation(),
      ! is not just a simple-minded time-mean of cos(sza) over the
      ! REFRESH interval, but a flux-weighted time-mean of cos(sza).
      mu0 = real(ZT, kind=wp)
      ! 2. total solar irradiance, tsi, NORMAL to solar beam [W/m2]
      ! NB: the naive approach would be to use a scalar, tsi = SC * DIST.
      ! But tsi is ultimately multiplied by mu0 in the radiative transfer,
      ! to yield a vertical flux, and then should be equal to the actual
      ! vertical solar flux at the TOA, namely SLR1D. We must therefore
      ! use SLR1D/mu0. This is not the same as SC * DIST, not only because
      ! of time co-variation of sza and DIST, but also because of the more
      ! nuanced definition of ZTH in item 1 above.
      ! Note: mu0 cannot be zero since running for daytime columns only.
      tsi = real(SLR1D, kind=wp) / mu0
      ! 3. surface albedos
      ! NIR bands (1-9: 820-12850 cm-1, 0.778-12.195 microns)
      do ib=1,9
        sfc_alb_dir(ib,:)  = real(ALBNR, kind=wp)
        sfc_alb_dif(ib,:)  = real(ALBNF, kind=wp)
      enddo
      ! UV/visible bands (11-14: 16000-50000 cm-1, 0.200-0.625 micron)
      do ib=11,14
        sfc_alb_dir(ib,:)  = real(ALBVR, kind=wp)
        sfc_alb_dif(ib,:)  = real(ALBVF, kind=wp)
      enddo
      ! Transition band (10, 12850-16000 cm-1, 0.625-0.778 micron)
      ! Take average, dmlee
      sfc_alb_dir(10,:) = real((ALBVR+ALBNR)/2., kind=wp)
      sfc_alb_dif(10,:) = real((ALBVF+ALBNF)/2., kind=wp)

      ! basic profiles
      p_lay = real(PL , kind=wp)
      t_lay = real(T  , kind=wp)
      p_lev = real(PLE, kind=wp)

      ! RRTMGP's rte_sw takes a vertical ordering flag
      ! (no need to flip columns as with RRTMG)
      top_at_1 = p_lay(1, 1) < p_lay(1, LM)
      _ASSERT(top_at_1, 'unexpected vertical ordering')

      ! layer pressure thicknesses used for cloud water path calculations
      ! (do before any KLUGE to top pressure so optical paths wont be affected)
      ! (also better to use these unKLUGED pressure intervals in t_lev calculation)
      dp_wp = p_lev(:,2:LM+1) - p_lev(:,1:LM)

      ! pmn: pressure KLUGE
      ! Because currently k_dist%press_ref_min ~ 1.005 > GEOS-5 ptop of 1.0 Pa.
      ! Find better solution, perhaps getting AER to add a higher top.
      press_ref_min = k_dist%get_press_min()
      ptop = minval(p_lev(:,1))
      if (press_ref_min > ptop) then
        ! allow a small increase of ptop
        if (press_ref_min - ptop <= ptop * ptop_increase_OK_fraction) then
          where (p_lev(:,1) < press_ref_min) p_lev(:,1) = press_ref_min
          ! make sure no pressure ordering issues were created
          _ASSERT(all(p_lev(:,1) < p_lay(:,1)), 'pressure kluge causes misordering')
        else
          write(*,*) ' A ', ptop_increase_OK_fraction, &
                       ' fractional increase of ptop was insufficient'
          write(*,*) ' RRTMGP, GEOS-5 top (Pa)', press_ref_min, ptop
          TEST_('Model top too high for RRTMGP')
        endif
      endif

      ! pmn: temperature KLUGE
      ! Currently k_dist%temp_ref_min = 160K but GEOS-5 has a global minimum
      ! temperature below this occasionally (< 1% of time). (The lowest temp
      ! seen so far is above 145K). Consequently we will limit min(t_lay) to
      ! 160K.
      ! Find better solution, perhaps getting AER to produce a table with a
      ! lower minimum temperature.
      temp_ref_min = k_dist%get_temp_min()
      tmin = minval(t_lay)
      if (temp_ref_min > tmin) then
        ! allow a small increase of tmin
        call MAPL_GetResource (MAPL, &
           tmin_increase_OK_Kelvin, 'RRTMGP_SW_TMIN_INC_OK_K:', &
           DEFAULT = 15._wp, __RC__)
        if (temp_ref_min - tmin <= tmin_increase_OK_Kelvin) then
          where (t_lay < temp_ref_min) t_lay = temp_ref_min
        else
          write(*,*) ' A ', tmin_increase_OK_Kelvin, &
                       'K increase of tmin was insufficient'
          write(*,*) ' RRTMGP, GEOS-5 t_min (K)', temp_ref_min, tmin
          TEST_('Found excessively cold model temperature for RRTMGP')
        endif
      endif

      ! dzmid(k) is separation [m] between midpoints of layers k and k+1 (sign not important, +ve here).
      ! dz ~ RT/g x dp/p by hydrostatic eqn and ideal gas eqn. The jump from LAYER k to k+1 is centered
      ! on LEVEL k+1 since the LEVEL indices are one-based.
      ! pmn: note that the dzmid calculation depends on the t_lev. Though t_lev is a temporary here, it
      ! is an important variable in the LW, where its calculation must occur after the t_lay KLUGE. So,
      ! for consistency with the LW, this t_lev and dzmid calculation is placed after the t_lay KLUGE.
      allocate(t_lev(ncol),__STAT__)
      do k = 1,LM-1
        ! t_lev are interior interface temperatures at level k+1
        t_lev = (t_lay(:,k) * dp_wp(:,k+1) + t_lay(:,k+1) * dp_wp(:,k)) / (dp_wp(:,k+1) + dp_wp(:,k))
        dzmid(:,k) = t_lev * real(MAPL_RGAS/MAPL_GRAV,kind=wp) * (p_lay(:,k+1) - p_lay(:,k)) / p_lev(:,k+1)
      end do
      deallocate(t_lev,__STAT__)

      ! allocation of output arrays
      allocate(flux_up_clrsky (ncol,LM+1), flux_net_clrsky(ncol,LM+1), __STAT__)
      allocate(flux_up_allsky (ncol,LM+1), flux_net_allsky(ncol,LM+1), __STAT__)
      allocate(bnd_flux_dn_allsky (ncol,LM+1,nbnd), &
               bnd_flux_net_allsky(ncol,LM+1,nbnd), &
               bnd_flux_dir_allsky(ncol,LM+1,nbnd), __STAT__)

      ! =====================================================================================
      ! IMPORTANT: Specify the type (#streams) of the SW RT calculations in optical_props
      ! =====================================================================================
      ! While the aerosol system currently provides two-stream properties, as do the cloud
      ! optics files, we may choose any number of streams for the actual RT calculations by
      ! the appropriate instantiation of optical_props here. The increment() statements below
      ! implicitly convert all component optical properties to this number of streams.
      ! Everything else in the code should adapt polymorphically without modification.
      ! Options are: 1scl (no scattering), 2str (2-stream), or nstr (n-stream).
      ! For nstr, must also specify the number of phase function moments (nmom) below.
      ! =====================================================================================

      ! instantiate optical_props with desired streams
      allocate(ty_optical_props_2str::optical_props,__STAT__)  ! <-- choose 2-stream SW

      ! initialize spectral discretiz'n and gpt mapping of optical_props
      TEST_(optical_props%init(k_dist))

      ! Used only if nstr (and then must be >= 2)
      nmom = 2

      ! get cloud optical properties (band-only)
      ! pmn: some of this could be done only once per run ...

      ! load and init cloud_optics from file:
      ! gets appropriate coefficients needed to calculate
      ! cloud optical properties from cloud physical properties
      call MAPL_GetResource( &
        MAPL, cloud_optics_file, "RRTMGP_CLOUD_OPTICS_SW:", &
        DEFAULT='rrtmgp-clouds-sw.nc', __RC__)
      call MAPL_GetResource( &
        MAPL, cloud_optics_type, "RRTMGP_CLOUD_OPTICS_TYPE_SW:", &
        DEFAULT='LUT', __RC__)
      call MAPL_TimerOn(MAPL,"--RRTMGP_IO_CLOUDS",__RC__)
      if (trim(cloud_optics_type)=='LUT') then
        call load_cld_lutcoeff (cloud_optics, cloud_optics_file)
      elseif (trim(cloud_optics_type)=='PADE') then
        call load_cld_padecoeff(cloud_optics, cloud_optics_file)
      else
        TEST_('unknown cloud_optics_type: '//trim(cloud_optics_file))
      end if
      call MAPL_TimerOff(MAPL,"--RRTMGP_IO_CLOUDS",__RC__)

      ! ice surface roughness category for Yang (2013) ice optics
      ! icergh: 1 = none, 2 = medium, 3 = high
      call MAPL_GetResource( &
        MAPL, icergh, "RRTMGP_ICE_ROUGHNESS_SW:", &
        DEFAULT=2, __RC__)
      TEST_(cloud_optics%set_ice_roughness(icergh))

      ! delta-scaling if of course applied by default
      ! ... you can turn it off for debugging purposes
      call MAPL_GetResource ( &
        MAPL, rrtmgp_delta_scale, LABEL='RRTMGP_DELTA_SCALE:', &
        DEFAULT=.TRUE., __RC__)
      call MAPL_GetResource ( &
        MAPL, rrtmgp_use_rrtmg_iceflg3_like_forwice, &
        LABEL='RRTMGP_USE_RRTMG_ICEFLG3_LIKE_FORWICE:', &
        DEFAULT=.TRUE., __RC__)

      ! cloud optics file is currently two-stream
      ! increment() will handle appropriate stream conversions
      allocate(ty_optical_props_2str::cloud_props_bnd_liq,__STAT__)
      allocate(ty_optical_props_2str::cloud_props_bnd_ice,__STAT__)

      ! band-only initialization for pre-mcICA cloud optical properties
      TEST_(cloud_props_bnd_liq%init(k_dist%get_band_lims_wavenumber()))
      TEST_(cloud_props_bnd_ice%init(k_dist%get_band_lims_wavenumber()))

      ! g-point version for McICA sampled cloud optical properties
      select type (cloud_props_bnd_liq)
        class is (ty_optical_props_2str)
          allocate(ty_optical_props_2str::cloud_props_gpt_liq,__STAT__)
        class default
          TEST_('cloud optical properties (liq) hardwired 2-stream for now')
      end select
      select type (cloud_props_bnd_ice)
        class is (ty_optical_props_2str)
          allocate(ty_optical_props_2str::cloud_props_gpt_ice,__STAT__)
        class default
          TEST_('cloud optical properties (ice) hardwired 2-stream for now')
      end select
      TEST_(cloud_props_gpt_liq%init(k_dist))
      TEST_(cloud_props_gpt_ice%init(k_dist))

      ! read desired cloud overlap type
      call MAPL_GetResource( &
        MAPL, cloud_overlap_type, "RRTMGP_CLOUD_OVERLAP_TYPE_SW:", &
        DEFAULT='GEN_MAX_RAN_OVERLAP', __RC__)

      ! GEN_MAX_RAN_OVERLAP uses correlation lengths
      !   and possibly inhomogeneous condensate
      gen_mro = (cloud_overlap_type == "GEN_MAX_RAN_OVERLAP")
      if (gen_mro) then

        ! condensate inhomogeneous?
        ! see RadiationGC initialization
        cond_inhomo = condensate_inhomogeneous()

        ! Compute decorrelation length scales [m]
        allocate(adl(ncol),__STAT__)
        call correlation_length_cloud_fraction(ncol, ncol, doy, alat, adl)
        if (cond_inhomo) then
          allocate(rdl(ncol),__STAT__)
          call correlation_length_condensate(ncol, ncol, doy, alat, rdl)
        endif

      endif

      ! ===============================================================================
      ! Random number setup:
      ! ===============================================================================
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
      !   To get a different set of random numbers for the LW, for example, either a
      ! key change or a counter advance will be needed.
      !
      ! Time Component of key:
      ! ~~~~~~~~~~~~~~~~~~~~~~
      ! 1. No need to update more frequently than once per SW refresh.
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
      ! ===============================================================================

      allocate(seeds(3),__STAT__) ! 2-word key plus word1 of counter

      ! seed(1), the column part (word1) of key is set later

      ! get time part (word2) of key
      call ESMF_TimeSet (ReferenceTime, yy=2000, mm=1, dd=1, __RC__)
      call ESMF_AlarmGet(ALARM, RINGINTERVAL=RefreshInterval, __RC__)
      seeds(2) = int((CurrTime - ReferenceTime) / RefreshInterval)

      ! for SW start at counter=65,536
      seeds(3) = 65536

      ! set up aerosol optical properties
      need_aer_optical_props = (include_aerosols .and. implements_aerosol_optics)
      if (need_aer_optical_props) then
        ! aerosol optics system is currently two-stream
        ! increment() will handle appropriate stream conversions
        allocate(ty_optical_props_2str::aer_props,__STAT__)
        ! band-only initialization
        TEST_(aer_props%init(k_dist%get_band_lims_wavenumber()))
      end if

      !-------------------------------------------------------!
      ! Loop over blocks of blockSize columns                 !
      !  - choose rrtmgp_blockSize for memory/time efficiency !
      !  - one possible partial block is done at the end      !
      !-------------------------------------------------------!

      call MAPL_GetResource( MAPL, &
        rrtmgp_blockSize, "RRTMGP_SW_BLOCKSIZE:", DEFAULT=4, __RC__)
      _ASSERT(rrtmgp_blockSize >= 1, 'bad RRTMGP_SW_BLOCKSIZE')

      ! for random numbers, for efficiency, reserve the maximum possible
      ! subset of columns (rrtmgp_blockSize) since column index is last
      allocate(urand(ngpt,LM,rrtmgp_blocksize),__STAT__)
      if (gen_mro) then
        allocate(urand_aux(ngpt,LM,rrtmgp_blocksize),__STAT__)
        if (cond_inhomo) then
          allocate(urand_cond    (ngpt,LM,rrtmgp_blocksize),__STAT__)
          allocate(urand_cond_aux(ngpt,LM,rrtmgp_blocksize),__STAT__)
        end if
      end if

      ! number of FULL blocks by integer division
      nBlocks = ncol/rrtmgp_blockSize

      ! allocate intermediate arrays for FULL blocks
      if (nBlocks > 0) then

        ! block size UNTIL possible final partial block
        ncols_block = rrtmgp_blockSize

        allocate(toa_flux(ncols_block,ngpt),__STAT__)
        allocate(cld_mask(ncols_block,LM,ngpt),__STAT__)
        allocate(forwliq(ncols_block,LM,ngpt),__STAT__)
        allocate(forwice(ncols_block,LM,ngpt),__STAT__)
        if (gen_mro) then
          allocate(alpha(ncols_block,LM-1),__STAT__)
          if (cond_inhomo) then
            allocate(rcorr(ncols_block,LM-1),__STAT__)
            allocate(zcw(ncols_block,LM,ngpt),__STAT__)
          endif
        endif
        if (include_aerosols) &
          allocate(ClearCounts(4,ncols_block),__STAT__)

        ! in-cloud cloud optical props
        select type (cloud_props_bnd_liq)
          class is (ty_optical_props_2str)
            TEST_(cloud_props_bnd_liq%alloc_2str(ncols_block,LM))
        end select
        select type (cloud_props_bnd_ice)
          class is (ty_optical_props_2str)
            TEST_(cloud_props_bnd_ice%alloc_2str(ncols_block,LM))
        end select
        select type (cloud_props_gpt_liq)
          class is (ty_optical_props_2str)
            TEST_(cloud_props_gpt_liq%alloc_2str(ncols_block,LM))
        end select
        select type (cloud_props_gpt_ice)
          class is (ty_optical_props_2str)
            TEST_(cloud_props_gpt_ice%alloc_2str(ncols_block,LM))
        end select

        ! aerosol optical props
        if (need_aer_optical_props) then
          select type (aer_props)
            class is (ty_optical_props_2str)
              TEST_(aer_props%alloc_2str(ncols_block,LM))
          end select
        end if

        ! gas+aer+cld optical properties
        select type (optical_props)
          class is (ty_optical_props_1scl)
            TEST_(optical_props%alloc_1scl(ncols_block,LM))
          class is (ty_optical_props_2str)
            TEST_(optical_props%alloc_2str(ncols_block,LM))
          class is (ty_optical_props_nstr)
            TEST_(optical_props%alloc_nstr(nmom,ncols_block,LM))
        end select

      end if

      ! add final partial block if necessary
      partial_block = mod(ncol,rrtmgp_blockSize) /= 0
      if (partial_block) then
        partial_blockSize = ncol - nBlocks * rrtmgp_blockSize
        nBlocks = nBlocks + 1
      endif

      ! loop over all blocks
      do b = 1,nBlocks

        ! only the FINAL block can be partial
        if (b == nBlocks .and. partial_block) then
          ncols_block = partial_blockSize

          if (b > 1) then
            ! one or more full blocks already processed
            deallocate(toa_flux,      __STAT__)
            deallocate(cld_mask,      __STAT__)
            deallocate(forwliq,       __STAT__)
            deallocate(forwice,       __STAT__)
            if (gen_mro) then
              deallocate(alpha,       __STAT__)
              if (cond_inhomo) then
                deallocate(rcorr,zcw, __STAT__)
              endif
            endif
            if (include_aerosols) &
              deallocate(ClearCounts, __STAT__)
          endif

          allocate(toa_flux(ncols_block,ngpt),    __STAT__)
          allocate(cld_mask(ncols_block,LM,ngpt), __STAT__)
          allocate(forwliq(ncols_block,LM,ngpt),  __STAT__)
          allocate(forwice(ncols_block,LM,ngpt),  __STAT__)
          if (gen_mro) then
            allocate(alpha(ncols_block,LM-1),     __STAT__)
            if (cond_inhomo) then
              allocate(rcorr(ncols_block,LM-1),   __STAT__)
              allocate(zcw(ncols_block,LM,ngpt),  __STAT__)
            endif
          endif
          if (include_aerosols) &
            allocate(ClearCounts(4,ncols_block),  __STAT__)

          ! ty_optical_props routines have an internal deallocation
          select type (cloud_props_bnd_liq)
            class is (ty_optical_props_2str)
              TEST_(cloud_props_bnd_liq%alloc_2str(ncols_block,LM))
          end select
          select type (cloud_props_bnd_ice)
            class is (ty_optical_props_2str)
              TEST_(cloud_props_bnd_ice%alloc_2str(ncols_block,LM))
          end select
          select type (cloud_props_gpt_liq)
            class is (ty_optical_props_2str)
              TEST_(cloud_props_gpt_liq%alloc_2str(ncols_block,LM))
          end select
          select type (cloud_props_gpt_ice)
            class is (ty_optical_props_2str)
              TEST_(cloud_props_gpt_ice%alloc_2str(ncols_block,LM))
          end select
          if (need_aer_optical_props) then
            select type (aer_props)
              class is (ty_optical_props_2str)
                TEST_(aer_props%alloc_2str(ncols_block,LM))
            end select
          end if
          select type (optical_props)
            class is (ty_optical_props_1scl)
              TEST_(optical_props%alloc_1scl(ncols_block,LM))
            class is (ty_optical_props_2str)
              TEST_(optical_props%alloc_2str(ncols_block,LM))
            class is (ty_optical_props_nstr)
              TEST_(optical_props%alloc_nstr(nmom,ncols_block,LM))
          end select

        endif  ! partial block

        ! prepare block
        colS = (b-1) * rrtmgp_blockSize + 1
        colE = colS + ncols_block - 1
        TEST_(gas_concs%get_subset(colS,ncols_block,gas_concs_block))

        call MAPL_TimerOn(MAPL,"--RRTMGP_GAS_OPTICS",__RC__)
        ! gas optics, including source functions
        error_msg = k_dist%gas_optics( &
          p_lay(colS:colE,:), p_lev(colS:colE,:), t_lay(colS:colE,:), &
          gas_concs_block, optical_props, toa_flux)
        TEST_(error_msg)
        call MAPL_TimerOff(MAPL,"--RRTMGP_GAS_OPTICS",__RC__)

        ! get block of aerosol optical props
        if (need_aer_optical_props) then
          select type (aer_props)
            class is (ty_optical_props_2str)

              ! load un-normalized optical properties from aerosol system
              aer_props%tau = real(TAUA(colS:colE,:,:),kind=wp)
              aer_props%ssa = real(SSAA(colS:colE,:,:),kind=wp)
              aer_props%g   = real(ASYA(colS:colE,:,:),kind=wp)

              ! renormalize
              where (aer_props%tau > 0._wp .and. aer_props%ssa > 0._wp)
                aer_props%g   = aer_props%g   / aer_props%ssa
                aer_props%ssa = aer_props%ssa / aer_props%tau
              elsewhere
                aer_props%tau = 0._wp
                aer_props%ssa = 0._wp
                aer_props%g   = 0._wp
              end where

            class default
              TEST_('aerosol optical properties hardwired 2-stream for now')
          end select
        end if

        call MAPL_TimerOn(MAPL,"--RRTMGP_CLOUD_OPTICS",__RC__)

        ! Make band in-cloud optical props from cloud_optics and mean in-cloud cloud water paths.
        ! These can be scaled later to account for sub-gridscale condensate inhomogeneity.
        ! Do phases separately to allow for different forward scattering, etc., per earlier note.
        ! liquid ...
        error_msg = cloud_optics%cloud_optics( &
          real(QQ3(colS:colE,:,2),kind=wp) * dp_wp(colS:colE,:) * cwp_fac, &  ! [g/m2]
          dummy_wp(colS:colE,:), & 
          min( max( real(RR3(colS:colE,:,2),kind=wp), &  ! [microns]
            cloud_optics%get_min_radius_liq()), &
            cloud_optics%get_max_radius_liq()), &
          dummy_wp(colS:colE,:), &
          cloud_props_bnd_liq)
        TEST_(error_msg)
        ! ice ...
        error_msg = cloud_optics%cloud_optics( &
          dummy_wp(colS:colE,:), &
          real(QQ3(colS:colE,:,1),kind=wp) * dp_wp(colS:colE,:) * cwp_fac, &  ! [g/m2]
          dummy_wp(colS:colE,:), &
          min( max( real(RR3(colS:colE,:,1),kind=wp), &  ! [microns]
            cloud_optics%get_min_radius_ice()), &
            cloud_optics%get_max_radius_ice()), &
          cloud_props_bnd_ice)
        TEST_(error_msg)

        call MAPL_TimerOff(MAPL,"--RRTMGP_CLOUD_OPTICS",__RC__)

        call MAPL_TimerOn(MAPL,"--RRTMGP_MCICA",__RC__)

!!TODO: need to resolve diff between prob of max vs ran and correlation coeff in both paper and code

        ! exponential inter-layer correlations
        ! [alpha|rcorr](k) is correlation between layers k and k+1
        ! dzmid(k) is separation between midpoints of layers k and k+1
        if (gen_mro) then
          do ilay = 1,LM-1
            ! cloud fraction correlation
            alpha(:,ilay) = exp(-abs(dzmid(colS:colE,ilay))/real(adl(colS:colE),kind=wp))
          enddo
          if (cond_inhomo) then
            do ilay = 1,LM-1
              ! condensate correlation
              rcorr(:,ilay) = exp(-abs(dzmid(colS:colE,ilay))/real(rdl(colS:colE),kind=wp))
            enddo
          endif
        endif

        ! generate McICA random numbers for block
        ! Perhaps later this can be parallelized?
#ifdef HAVE_MKL
        do isub = 1, ncols_block
          ! local 1d column index
          icol = colS + isub - 1
          ! initialize the Philox PRNG
          ! set word1 of key based on GLOBAL location
          ! 32-bits can hold all forseeable resolutions
          seeds(1) = nint(Jg1D(icol)) * IM_World + nint(Ig1D(icol))
          ! instantiate a random number stream for the column
          call rng%init(VSL_BRNG_PHILOX4X32X10,seeds)
          ! draw the random numbers for the column
          urand(:,:,isub) = reshape(rng%get_random(ngpt*LM),(/ngpt,LM/))
          if (gen_mro) then
            urand_aux(:,:,isub) = reshape(rng%get_random(ngpt*LM),(/ngpt,LM/))
            if (cond_inhomo) then
              urand_cond    (:,:,isub) = reshape(rng%get_random(ngpt*LM),(/ngpt,LM/))
              urand_cond_aux(:,:,isub) = reshape(rng%get_random(ngpt*LM),(/ngpt,LM/))
            endif
          end if
          ! free the rng
          call rng%end()
        end do
#endif

        ! cloud sampling to gpoints
        select case (cloud_overlap_type)
          case ("MAX_RAN_OVERLAP")
            error_msg = sampled_mask_max_ran( &
              urand(:,:,1:ncols_block), real(CL(colS:colE,:),kind=wp), cld_mask)
            TEST_(error_msg)
          case ("EXP_RAN_OVERLAP")
            ! corr_coeff(ncols_block,LM-1) is an inter-layer correlation coefficient
            ! to be provided ... it is not the same as alpha, which is a probability
!           error_msg = sampled_mask_exp_ran( &
!             urand(:,:,1:ncols_block), real(CL(colS:colE,:),kind=wp), corr_coeff, cld_mask)
!           TEST_(error_msg)
            TEST_('EXP_RAN_OVERLAP not implemented yet')
          case ("GEN_MAX_RAN_OVERLAP")
            ! a scheme like Oreopoulos et al. 2012 (doi:10.5194/acp-12-9097-2012) in which both
            ! cloud presence and cloud condensate are separately generalized maximum-random:
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
                cld_frac = CL(icol,ilay)

                ! if grid-box clear, no subgrid variability
                if (cld_frac <= 0.) then
                  cld_mask(isub,ilay,:) = .false.
                else
                  ! subgrid-scale cloud mask
                  cld_mask(isub,ilay,:) = urand(:,ilay,isub) < cld_frac

                  ! subgrid-scale condensate
                  if (cond_inhomo) then
                    ! level of condensate inhomogeneity based on cloud fraction.
                    if (cld_frac > 0.99) then
                      sigma_qcw = 0.5
                    elseif (cld_frac > 0.9) then
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
            TEST_('RRTMGP_SW: unknown cloud overlap')
        end select

        ! draw McICA optical property samples (band->gpt)
        TEST_(draw_samples(cld_mask, cloud_props_bnd_liq, cloud_props_gpt_liq))
        TEST_(draw_samples(cld_mask, cloud_props_bnd_ice, cloud_props_gpt_ice))

        ! Scaling to sub-gridscale water paths:
        ! since tau for each phase is linear in the phase's water path
        ! and since the scaling zcw applies equally to both phases, the
        ! total g-point optical thickness tau will scale with zcw.
        if (gen_mro) then
          if (cond_inhomo) &
            where (cld_mask) cloud_props_gpt_liq%tau = cloud_props_gpt_liq%tau * zcw
            where (cld_mask) cloud_props_gpt_ice%tau = cloud_props_gpt_ice%tau * zcw
        end if

        call MAPL_TimerOff(MAPL,"--RRTMGP_MCICA",__RC__)

        ! REFRESH super-layer diagnostics (before delta-scaling TAUs).
        ! ** Calculated from subcolumn ensemble, so stochastic **
        ! -------------------------------------------------------
        call MAPL_TimerOn(MAPL,"--RRTMGP_SPRLYR_DIAGS",__RC__)
        if (include_aerosols) then

          ! super-layer cloud fractions
          call clearCounts_threeBand( &
            ncols_block, ncols_block, ngpt, LM, LCLDLM, LCLDMH, &
            reshape(cld_mask,[LM,ngpt,ncols_block],order=[3,1,2]), &
            ClearCounts)
          do isub = 1,ncols_block
            icol = colS + isub - 1
            CLDTS(icol) = 1. - ClearCounts(1,isub)/float(ngpt)
            CLDHS(icol) = 1. - ClearCounts(2,isub)/float(ngpt)
            CLDMS(icol) = 1. - ClearCounts(3,isub)/float(ngpt)
            CLDLS(icol) = 1. - ClearCounts(4,isub)/float(ngpt)
          end do

          ! in-cloud optical thicknesses in PAR super-band
          ! (weighted across and within bands by TOA incident flux)
          do isub = 1,ncols_block
            icol = colS + isub - 1

#ifdef SOLAR_RADVAL
            ! default (no cloud) for TAUx variant 
            TAUTP(icol) = 0.
            TAUHP(icol) = 0.
            TAUMP(icol) = 0.
            TAULP(icol) = 0.
#endif

            ! default (no cloud) for COTx variant
            COTTP(icol) = MAPL_UNDEF
            COTHP(icol) = MAPL_UNDEF
            COTMP(icol) = MAPL_UNDEF
            COTLP(icol) = MAPL_UNDEF

            ! zero denom- and numerator accumulators
            COTDTP(icol) = 0.; COTNTP(icol) = 0.
            COTDHP(icol) = 0.; COTNHP(icol) = 0.
            COTDMP(icol) = 0.; COTNMP(icol) = 0.
            COTDLP(icol) = 0.; COTNLP(icol) = 0.
#ifdef SOLAR_RADVAL
            COTLDTP(icol) = 0.; COTLNTP(icol) = 0.; COTIDTP(icol) = 0.; COTINTP(icol) = 0.
            COTLDHP(icol) = 0.; COTLNHP(icol) = 0.; COTIDHP(icol) = 0.; COTINHP(icol) = 0.
            COTLDMP(icol) = 0.; COTLNMP(icol) = 0.; COTIDMP(icol) = 0.; COTINMP(icol) = 0.
            COTLDLP(icol) = 0.; COTLNLP(icol) = 0.; COTIDLP(icol) = 0.; COTINLP(icol) = 0.
            SSALDTP(icol) = 0.; SSALNTP(icol) = 0.; SSAIDTP(icol) = 0.; SSAINTP(icol) = 0.
            SSALDHP(icol) = 0.; SSALNHP(icol) = 0.; SSAIDHP(icol) = 0.; SSAINHP(icol) = 0.
            SSALDMP(icol) = 0.; SSALNMP(icol) = 0.; SSAIDMP(icol) = 0.; SSAINMP(icol) = 0.
            SSALDLP(icol) = 0.; SSALNLP(icol) = 0.; SSAIDLP(icol) = 0.; SSAINLP(icol) = 0.
            ASMLDTP(icol) = 0.; ASMLNTP(icol) = 0.; ASMIDTP(icol) = 0.; ASMINTP(icol) = 0.
            ASMLDHP(icol) = 0.; ASMLNHP(icol) = 0.; ASMIDHP(icol) = 0.; ASMINHP(icol) = 0.
            ASMLDMP(icol) = 0.; ASMLNMP(icol) = 0.; ASMIDMP(icol) = 0.; ASMINMP(icol) = 0.
            ASMLDLP(icol) = 0.; ASMLNLP(icol) = 0.; ASMIDLP(icol) = 0.; ASMINLP(icol) = 0.
#endif

            ! can only be non-zero for potentially cloudy columns
            if (any(CL(icol,:) > 0.)) then

              ! accumulate over gpts/subcolumns
              do ib = 1, nbnd
                do igpt = band_lims_gpt(1,ib), band_lims_gpt(2,ib)
       
                  ! band weights for photosynthetically active radiation (PAR)
                  ! Bands 11-12 (0.345-0.625 um) plus half transition band 10 (0.625-0.778 um)
                  if (ib >= 11 .and. ib <= 12) then
                    wgt = 1.0
                  else if (ib == 10) then
                    wgt = 0.5
                  else
                    ! no contribution to PAR
                    cycle
                  end if

                  ! TOA flux weighting
                  ! (note: neither the adjustment of toa_flux to our tsi
                  ! or for zenith angle are needed yet since this weighting
                  ! is over gpoint and is normalized for EACH icol)
                  wgt = wgt * toa_flux(isub,igpt)

                  ! low pressure layer
                  sltaulp = sum(cloud_props_gpt_liq%tau(isub,LCLDLM:LM,igpt))
                  sitaulp = sum(cloud_props_gpt_ice%tau(isub,LCLDLM:LM,igpt))
                  staulp = sltaulp + sitaulp
                  if (staulp > 0.) then
                    COTDLP(icol) = COTDLP(icol) + wgt
                    COTNLP(icol) = COTNLP(icol) + wgt * staulp
                  end if
#ifdef SOLAR_RADVAL
                  sltaussalp = 0.; sltaussaglp = 0.
                  if (sltaulp > 0.) then
                    select type(cloud_props_gpt_liq)
                    class is (ty_optical_props_2str)
                      sltaussalp = sum(cloud_props_gpt_liq%tau(isub,LCLDLM:LM,igpt) * &
                                       cloud_props_gpt_liq%ssa(isub,LCLDLM:LM,igpt))
                      sltaussaglp = sum(cloud_props_gpt_liq%tau(isub,LCLDLM:LM,igpt) * &
                                        cloud_props_gpt_liq%ssa(isub,LCLDLM:LM,igpt) * &
                                        cloud_props_gpt_liq%g  (isub,LCLDLM:LM,igpt))
                    end select
                    COTLDLP(icol) = COTLDLP(icol) + wgt
                    COTLNLP(icol) = COTLNLP(icol) + wgt * sltaulp
                    SSALDLP(icol) = SSALDLP(icol) + wgt * sltaulp
                    SSALNLP(icol) = SSALNLP(icol) + wgt * sltaussalp
                    ASMLDLP(icol) = ASMLDLP(icol) + wgt * sltaussalp
                    ASMLNLP(icol) = ASMLNLP(icol) + wgt * sltaussaglp
                  end if
                  sitaussalp = 0.; sitaussaglp = 0.
                  if (sitaulp > 0.) then
                    select type(cloud_props_gpt_ice)
                    class is (ty_optical_props_2str)
                      sitaussalp = sum(cloud_props_gpt_ice%tau(isub,LCLDLM:LM,igpt) * &
                                       cloud_props_gpt_ice%ssa(isub,LCLDLM:LM,igpt))
                      sitaussaglp = sum(cloud_props_gpt_ice%tau(isub,LCLDLM:LM,igpt) * &
                                        cloud_props_gpt_ice%ssa(isub,LCLDLM:LM,igpt) * &
                                        cloud_props_gpt_ice%g  (isub,LCLDLM:LM,igpt))
                    end select
                    COTIDLP(icol) = COTIDLP(icol) + wgt
                    COTINLP(icol) = COTINLP(icol) + wgt * sitaulp
                    SSAIDLP(icol) = SSAIDLP(icol) + wgt * sitaulp
                    SSAINLP(icol) = SSAINLP(icol) + wgt * sitaussalp
                    ASMIDLP(icol) = ASMIDLP(icol) + wgt * sitaussalp
                    ASMINLP(icol) = ASMINLP(icol) + wgt * sitaussaglp
                  end if
#endif

                  ! mid pressure layer
                  sltaump = sum(cloud_props_gpt_liq%tau(isub,LCLDMH:LCLDLM-1,igpt))
                  sitaump = sum(cloud_props_gpt_ice%tau(isub,LCLDMH:LCLDLM-1,igpt))
                  staump = sltaump + sitaump
                  if (staump > 0.) then
                    COTDMP(icol) = COTDMP(icol) + wgt
                    COTNMP(icol) = COTNMP(icol) + wgt * staump
                  end if
#ifdef SOLAR_RADVAL
                  sltaussamp = 0.; sltaussagmp = 0.
                  if (sltaump > 0.) then
                    select type(cloud_props_gpt_liq)
                    class is (ty_optical_props_2str)
                      sltaussamp = sum(cloud_props_gpt_liq%tau(isub,LCLDMH:LCLDLM-1,igpt) * &
                                       cloud_props_gpt_liq%ssa(isub,LCLDMH:LCLDLM-1,igpt))
                      sltaussagmp = sum(cloud_props_gpt_liq%tau(isub,LCLDMH:LCLDLM-1,igpt) * &
                                        cloud_props_gpt_liq%ssa(isub,LCLDMH:LCLDLM-1,igpt) * &
                                        cloud_props_gpt_liq%g  (isub,LCLDMH:LCLDLM-1,igpt))
                    end select
                    COTLDMP(icol) = COTLDMP(icol) + wgt
                    COTLNMP(icol) = COTLNMP(icol) + wgt * sltaump
                    SSALDMP(icol) = SSALDMP(icol) + wgt * sltaump
                    SSALNMP(icol) = SSALNMP(icol) + wgt * sltaussamp
                    ASMLDMP(icol) = ASMLDMP(icol) + wgt * sltaussamp
                    ASMLNMP(icol) = ASMLNMP(icol) + wgt * sltaussagmp
                  end if
                  sitaussamp = 0.; sitaussagmp = 0.
                  if (sitaump > 0.) then
                    select type(cloud_props_gpt_ice)
                    class is (ty_optical_props_2str)
                      sitaussamp = sum(cloud_props_gpt_ice%tau(isub,LCLDMH:LCLDLM-1,igpt) * &
                                       cloud_props_gpt_ice%ssa(isub,LCLDMH:LCLDLM-1,igpt))
                      sitaussagmp = sum(cloud_props_gpt_ice%tau(isub,LCLDMH:LCLDLM-1,igpt) * &
                                        cloud_props_gpt_ice%ssa(isub,LCLDMH:LCLDLM-1,igpt) * &
                                        cloud_props_gpt_ice%g  (isub,LCLDMH:LCLDLM-1,igpt))
                    end select
                    COTIDMP(icol) = COTIDMP(icol) + wgt
                    COTINMP(icol) = COTINMP(icol) + wgt * sitaump
                    SSAIDMP(icol) = SSAIDMP(icol) + wgt * sitaump
                    SSAINMP(icol) = SSAINMP(icol) + wgt * sitaussamp
                    ASMIDMP(icol) = ASMIDMP(icol) + wgt * sitaussamp
                    ASMINMP(icol) = ASMINMP(icol) + wgt * sitaussagmp
                  end if
#endif

                  ! high pressure layer
                  sltauhp = sum(cloud_props_gpt_liq%tau(isub,1:LCLDMH-1,igpt))
                  sitauhp = sum(cloud_props_gpt_ice%tau(isub,1:LCLDMH-1,igpt))
                  stauhp = sltauhp + sitauhp
                  if (stauhp > 0.) then
                    COTDHP(icol) = COTDHP(icol) + wgt
                    COTNHP(icol) = COTNHP(icol) + wgt * stauhp
                  end if
#ifdef SOLAR_RADVAL
                  sltaussahp = 0.; sltaussaghp = 0.
                  if (sltauhp > 0.) then
                    select type(cloud_props_gpt_liq)
                    class is (ty_optical_props_2str)
                      sltaussahp = sum(cloud_props_gpt_liq%tau(isub,1:LCLDMH-1,igpt) * &
                                       cloud_props_gpt_liq%ssa(isub,1:LCLDMH-1,igpt))
                      sltaussaghp = sum(cloud_props_gpt_liq%tau(isub,1:LCLDMH-1,igpt) * &
                                        cloud_props_gpt_liq%ssa(isub,1:LCLDMH-1,igpt) * &
                                        cloud_props_gpt_liq%g  (isub,1:LCLDMH-1,igpt))
                    end select
                    COTLDHP(icol) = COTLDHP(icol) + wgt
                    COTLNHP(icol) = COTLNHP(icol) + wgt * sltauhp
                    SSALDHP(icol) = SSALDHP(icol) + wgt * sltauhp
                    SSALNHP(icol) = SSALNHP(icol) + wgt * sltaussahp
                    ASMLDHP(icol) = ASMLDHP(icol) + wgt * sltaussahp
                    ASMLNHP(icol) = ASMLNHP(icol) + wgt * sltaussaghp
                  end if
                  sitaussahp = 0.; sitaussaghp = 0.
                  if (sitauhp > 0.) then
                    select type(cloud_props_gpt_ice)
                    class is (ty_optical_props_2str)
                      sitaussahp = sum(cloud_props_gpt_ice%tau(isub,1:LCLDMH-1,igpt) * &
                                       cloud_props_gpt_ice%ssa(isub,1:LCLDMH-1,igpt))
                      sitaussaghp = sum(cloud_props_gpt_ice%tau(isub,1:LCLDMH-1,igpt) * &
                                        cloud_props_gpt_ice%ssa(isub,1:LCLDMH-1,igpt) * &
                                        cloud_props_gpt_ice%g  (isub,1:LCLDMH-1,igpt))
                    end select
                    COTIDHP(icol) = COTIDHP(icol) + wgt
                    COTINHP(icol) = COTINHP(icol) + wgt * sitauhp
                    SSAIDHP(icol) = SSAIDHP(icol) + wgt * sitauhp
                    SSAINHP(icol) = SSAINHP(icol) + wgt * sitaussahp
                    ASMIDHP(icol) = ASMIDHP(icol) + wgt * sitaussahp
                    ASMINHP(icol) = ASMINHP(icol) + wgt * sitaussaghp
                  end if
#endif

                  ! whole subcolumn
                  sltautp = sltaulp + sltaump + sltauhp
                  sitautp = sitaulp + sitaump + sitauhp
                  stautp = staulp + staump + stauhp
                  if (stautp > 0.) then
                    COTDTP(icol) = COTDTP(icol) + wgt
                    COTNTP(icol) = COTNTP(icol) + wgt * stautp
                  end if
#ifdef SOLAR_RADVAL
                  sltaussatp = sltaussalp + sltaussamp + sltaussahp
                  sltaussagtp = sltaussaglp + sltaussagmp + sltaussaghp
                  if (sltautp > 0.) then
                    COTLDTP(icol) = COTLDTP(icol) + wgt
                    COTLNTP(icol) = COTLNTP(icol) + wgt * sltautp
                    SSALDTP(icol) = SSALDTP(icol) + wgt * sltautp
                    SSALNTP(icol) = SSALNTP(icol) + wgt * sltaussatp
                    ASMLDTP(icol) = ASMLDTP(icol) + wgt * sltaussatp
                    ASMLNTP(icol) = ASMLNTP(icol) + wgt * sltaussagtp
                  end if
                  sitaussatp = sitaussalp + sitaussamp + sitaussahp
                  sitaussagtp = sitaussaglp + sitaussagmp + sitaussaghp
                  if (sitautp > 0.) then
                    COTIDTP(icol) = COTIDTP(icol) + wgt
                    COTINTP(icol) = COTINTP(icol) + wgt * sitautp
                    SSAIDTP(icol) = SSAIDTP(icol) + wgt * sitautp
                    SSAINTP(icol) = SSAINTP(icol) + wgt * sitaussatp
                    ASMIDTP(icol) = ASMIDTP(icol) + wgt * sitaussatp
                    ASMINTP(icol) = ASMINTP(icol) + wgt * sitaussagtp
                  end if
#endif

                end do ! igpt
              end do ! ib

              ! normalize
              ! Note: TAUx defaults zero, COTx defaults MAPL_UNDEF
              if (COTDTP(icol) > 0. .and. COTNTP(icol) > 0.) then
                COTTP(icol) = COTNTP(icol) / COTDTP(icol)
#ifdef SOLAR_RADVAL
                TAUTP(icol) = COTTP(icol)
#endif
              end if

              if (COTDHP(icol) > 0. .and. COTNHP(icol) > 0.) then
                COTHP(icol) = COTNHP(icol) / COTDHP(icol)
#ifdef SOLAR_RADVAL
                TAUHP(icol) = COTHP(icol)
#endif
              end if

              if (COTDMP(icol) > 0. .and. COTNMP(icol) > 0.) then
                COTMP(icol) = COTNMP(icol) / COTDMP(icol)
#ifdef SOLAR_RADVAL
                TAUMP(icol) = COTMP(icol)
#endif
              end if

              if (COTDLP(icol) > 0. .and. COTNLP(icol) > 0.) then
                COTLP(icol) = COTNLP(icol) / COTDLP(icol)
#ifdef SOLAR_RADVAL
                TAULP(icol) = COTLP(icol)
#endif
              end if

            end if  ! potentially cloudy column 
          end do  ! isub
        end if  ! include_aerosols
        call MAPL_TimerOff(MAPL,"--RRTMGP_SPRLYR_DIAGS",__RC__)

        ! delta-scaling of cloud optical properties (accounts for forward scattering)
        call MAPL_TimerOn(MAPL,"--RRTMGP_DELTA_SCALE",__RC__)
        forwliq = 0.; forwice = 0.  ! default for no delta-scaling
        if (rrtmgp_delta_scale) then

          ! default delta-scaling for liquid
          select type(cloud_props_gpt_liq)
          class is (ty_optical_props_2str)
            forwliq = cloud_props_gpt_liq%g ** 2
          end select
          TEST_(cloud_props_gpt_liq%delta_scale(forwliq))

          if (rrtmgp_use_rrtmg_iceflg3_like_forwice) then
            ! non-default delta-scaling for ice (as in RRTMG iceflag==3)
            select type(cloud_props_gpt_ice)
            class is (ty_optical_props_2str)
              radice_lwr = cloud_optics%get_min_radius_ice()
              radice_upr = cloud_optics%get_max_radius_ice()
              do isub = 1,ncols_block
                icol = colS + isub - 1
                do ilay = 1,LM
                  ! only if at least potentially cloudy ...
                  if (CL(icol,ilay) > 0.) then
  
                    ! prepare for radice interpolation ...
                    ! first get radice consistent with RRTMGP ice cloud optics
                    radice = min(max(real(RR3(icol,ilay,1),kind=wp),radice_lwr),radice_upr)
                    ! now force into RRTMG's iceflag==3 reice binning range [5,140]um.
                    radice = min(max(radice,5._wp),140._wp)
                    ! RRTMG has 46 reice bins with 5um->radidx==1, 140um->radidx==46,
                    ! but radidx is forced to [1,45] so LIN2_ARG1 interpolation works.
                    radfac = (radice - 2._wp) / 3._wp
                    radidx = min(max(int(radfac),1),45)
                    rfint = radfac - real(radidx,kind=wp)
  
                    do ib = 1,nbnd
                      ! interpolate fdelta in radice for band ib
                      fdelta = LIN2_ARG1(fdlice3_rrtmgp,radidx,ib,rfint)
  
                      ! forwice calc for each g-point
                      do igpt = band_lims_gpt(1,ib),band_lims_gpt(2,ib)
                        if (cloud_props_gpt_ice%tau(isub,ilay,igpt) > 0.) then
                          forwice(isub,ilay,igpt) = min( &
                             fdelta + 0.5_wp / cloud_props_gpt_ice%ssa(isub,ilay,igpt), &
                             cloud_props_gpt_ice%g(isub,ilay,igpt))
                        endif
                      enddo  ! g-points
                    enddo  ! bands
  
                  endif  ! potentially cloudy
                enddo  ! layers
              enddo  ! columns
            end select
            TEST_(cloud_props_gpt_ice%delta_scale(forwice))
          else
            ! default delta-scaling for ice
            select type(cloud_props_gpt_ice)
            class is (ty_optical_props_2str)
              forwice = cloud_props_gpt_ice%g ** 2
            end select
            TEST_(cloud_props_gpt_ice%delta_scale(forwice))
          endif
        endif
        call MAPL_TimerOff(MAPL,"--RRTMGP_DELTA_SCALE",__RC__)

#ifdef SOLAR_RADVAL
        ! REFRESH super-layer diagnostics (after delta-scaling TAUs).
        ! ** Calculated from subcolumn ensemble, so stochastic **
        ! -------------------------------------------------------
        call MAPL_TimerOn(MAPL,"--RRTMGP_SPRLYR_DIAGS",__RC__)
        if (include_aerosols) then

          ! in-cloud optical thicknesses in PAR super-band
          ! (weighted across and within bands by TOA incident flux)
          do isub = 1,ncols_block
            icol = colS + isub - 1

            ! zero denom- and numerator accumulators
            CDSDTP(icol) = 0.; CDSNTP(icol) = 0.
            CDSDHP(icol) = 0.; CDSNHP(icol) = 0.
            CDSDMP(icol) = 0.; CDSNMP(icol) = 0.
            CDSDLP(icol) = 0.; CDSNLP(icol) = 0.

            CDSLDTP(icol) = 0.; CDSLNTP(icol) = 0.; CDSIDTP(icol) = 0.; CDSINTP(icol) = 0.
            CDSLDHP(icol) = 0.; CDSLNHP(icol) = 0.; CDSIDHP(icol) = 0.; CDSINHP(icol) = 0.
            CDSLDMP(icol) = 0.; CDSLNMP(icol) = 0.; CDSIDMP(icol) = 0.; CDSINMP(icol) = 0.
            CDSLDLP(icol) = 0.; CDSLNLP(icol) = 0.; CDSIDLP(icol) = 0.; CDSINLP(icol) = 0.

            SDSLDTP(icol) = 0.; SDSLNTP(icol) = 0.; SDSIDTP(icol) = 0.; SDSINTP(icol) = 0.
            SDSLDHP(icol) = 0.; SDSLNHP(icol) = 0.; SDSIDHP(icol) = 0.; SDSINHP(icol) = 0.
            SDSLDMP(icol) = 0.; SDSLNMP(icol) = 0.; SDSIDMP(icol) = 0.; SDSINMP(icol) = 0.
            SDSLDLP(icol) = 0.; SDSLNLP(icol) = 0.; SDSIDLP(icol) = 0.; SDSINLP(icol) = 0.

            ADSLDTP(icol) = 0.; ADSLNTP(icol) = 0.; ADSIDTP(icol) = 0.; ADSINTP(icol) = 0.
            ADSLDHP(icol) = 0.; ADSLNHP(icol) = 0.; ADSIDHP(icol) = 0.; ADSINHP(icol) = 0.
            ADSLDMP(icol) = 0.; ADSLNMP(icol) = 0.; ADSIDMP(icol) = 0.; ADSINMP(icol) = 0.
            ADSLDLP(icol) = 0.; ADSLNLP(icol) = 0.; ADSIDLP(icol) = 0.; ADSINLP(icol) = 0.

            FORLDTP(icol) = 0.; FORLNTP(icol) = 0.; FORIDTP(icol) = 0.; FORINTP(icol) = 0.
            FORLDHP(icol) = 0.; FORLNHP(icol) = 0.; FORIDHP(icol) = 0.; FORINHP(icol) = 0.
            FORLDMP(icol) = 0.; FORLNMP(icol) = 0.; FORIDMP(icol) = 0.; FORINMP(icol) = 0.
            FORLDLP(icol) = 0.; FORLNLP(icol) = 0.; FORIDLP(icol) = 0.; FORINLP(icol) = 0.

            ! can only be non-zero for potentially cloudy columns
            if (any(CL(icol,:) > 0.)) then

              ! accumulate over gpts/subcolumns
              do ib = 1, nbnd
                do igpt = band_lims_gpt(1,ib), band_lims_gpt(2,ib)

                  ! band weights for photosynthetically active radiation (PAR)
                  ! Bands 11-12 (0.345-0.625 um) plus half transition band 10 (0.625-0.778 um)
                  if (ib >= 11 .and. ib <= 12) then
                    wgt = 1.0
                  else if (ib == 10) then
                    wgt = 0.5
                  else
                    ! no contribution to PAR
                    cycle
                  end if

                  ! TOA flux weighting
                  ! (note: neither the adjustment of toa_flux to our tsi
                  ! or for zenith angle are needed yet since this weighting
                  ! is over gpoint and is normalized for EACH icol)
                  wgt = wgt * toa_flux(isub,igpt)

                  ! low pressure layer
                  sltaulp = sum(cloud_props_gpt_liq%tau(isub,LCLDLM:LM,igpt))
                  sitaulp = sum(cloud_props_gpt_ice%tau(isub,LCLDLM:LM,igpt))
                  staulp = sltaulp + sitaulp
                  if (staulp > 0.) then
                    CDSNLP(icol) = CDSNLP(icol) + wgt * staulp
                    CDSDLP(icol) = CDSDLP(icol) + wgt
                  end if
                  sltaussalp = 0.; sltaussaglp = 0.; sltaussaflp = 0.
                  if (sltaulp > 0.) then
                    select type(cloud_props_gpt_liq)
                    class is (ty_optical_props_2str)
                      sltaussalp = sum(cloud_props_gpt_liq%tau(isub,LCLDLM:LM,igpt) * &
                                       cloud_props_gpt_liq%ssa(isub,LCLDLM:LM,igpt))
                      sltaussaglp = sum(cloud_props_gpt_liq%tau(isub,LCLDLM:LM,igpt) * &
                                        cloud_props_gpt_liq%ssa(isub,LCLDLM:LM,igpt) * &
                                        cloud_props_gpt_liq%g  (isub,LCLDLM:LM,igpt))
                      sltaussaflp = sum(cloud_props_gpt_liq%tau(isub,LCLDLM:LM,igpt) * &
                                        cloud_props_gpt_liq%ssa(isub,LCLDLM:LM,igpt) * &
                                                        forwliq(isub,LCLDLM:LM,igpt))
                    end select
                    CDSLDLP(icol) = CDSLDLP(icol) + wgt
                    CDSLNLP(icol) = CDSLNLP(icol) + wgt * sltaulp
                    SDSLDLP(icol) = SDSLDLP(icol) + wgt * sltaulp
                    SDSLNLP(icol) = SDSLNLP(icol) + wgt * sltaussalp
                    ADSLDLP(icol) = ADSLDLP(icol) + wgt * sltaussalp
                    ADSLNLP(icol) = ADSLNLP(icol) + wgt * sltaussaglp
                    FORLDLP(icol) = FORLDLP(icol) + wgt * sltaussalp
                    FORLNLP(icol) = FORLNLP(icol) + wgt * sltaussaflp
                  end if
                  sitaussalp = 0.; sitaussaglp = 0.; sitaussaflp = 0.
                  if (sitaulp > 0.) then
                    select type(cloud_props_gpt_ice)
                    class is (ty_optical_props_2str)
                      sitaussalp = sum(cloud_props_gpt_ice%tau(isub,LCLDLM:LM,igpt) * &
                                       cloud_props_gpt_ice%ssa(isub,LCLDLM:LM,igpt))
                      sitaussaglp = sum(cloud_props_gpt_ice%tau(isub,LCLDLM:LM,igpt) * &
                                        cloud_props_gpt_ice%ssa(isub,LCLDLM:LM,igpt) * &
                                        cloud_props_gpt_ice%g  (isub,LCLDLM:LM,igpt))
                      sitaussaflp = sum(cloud_props_gpt_ice%tau(isub,LCLDLM:LM,igpt) * &
                                        cloud_props_gpt_ice%ssa(isub,LCLDLM:LM,igpt) * &
                                                        forwice(isub,LCLDLM:LM,igpt))
                    end select
                    CDSIDLP(icol) = CDSIDLP(icol) + wgt
                    CDSINLP(icol) = CDSINLP(icol) + wgt * sitaulp
                    SDSIDLP(icol) = SDSIDLP(icol) + wgt * sitaulp
                    SDSINLP(icol) = SDSINLP(icol) + wgt * sitaussalp
                    ADSIDLP(icol) = ADSIDLP(icol) + wgt * sitaussalp
                    ADSINLP(icol) = ADSINLP(icol) + wgt * sitaussaglp
                    FORIDLP(icol) = FORIDLP(icol) + wgt * sitaussalp
                    FORINLP(icol) = FORINLP(icol) + wgt * sitaussaflp
                  end if

                  ! mid pressure layer
                  sltaump = sum(cloud_props_gpt_liq%tau(isub,LCLDMH:LCLDLM-1,igpt))
                  sitaump = sum(cloud_props_gpt_ice%tau(isub,LCLDMH:LCLDLM-1,igpt))
                  staump = sltaump + sitaump
                  if (staump > 0.) then
                    CDSNMP(icol) = CDSNMP(icol) + wgt * staump
                    CDSDMP(icol) = CDSDMP(icol) + wgt
                  end if
                  sltaussamp = 0.; sltaussagmp = 0.; sltaussafmp = 0.
                  if (sltaump > 0.) then
                    select type(cloud_props_gpt_liq)
                    class is (ty_optical_props_2str)
                      sltaussamp = sum(cloud_props_gpt_liq%tau(isub,LCLDMH:LCLDLM-1,igpt) * &
                                       cloud_props_gpt_liq%ssa(isub,LCLDMH:LCLDLM-1,igpt))
                      sltaussagmp = sum(cloud_props_gpt_liq%tau(isub,LCLDMH:LCLDLM-1,igpt) * &
                                        cloud_props_gpt_liq%ssa(isub,LCLDMH:LCLDLM-1,igpt) * &
                                        cloud_props_gpt_liq%g  (isub,LCLDMH:LCLDLM-1,igpt))
                      sltaussafmp = sum(cloud_props_gpt_liq%tau(isub,LCLDMH:LCLDLM-1,igpt) * &
                                        cloud_props_gpt_liq%ssa(isub,LCLDMH:LCLDLM-1,igpt) * &
                                                        forwliq(isub,LCLDMH:LCLDLM-1,igpt))
                    end select
                    CDSLDMP(icol) = CDSLDMP(icol) + wgt
                    CDSLNMP(icol) = CDSLNMP(icol) + wgt * sltaump
                    SDSLDMP(icol) = SDSLDMP(icol) + wgt * sltaump
                    SDSLNMP(icol) = SDSLNMP(icol) + wgt * sltaussamp
                    ADSLDMP(icol) = ADSLDMP(icol) + wgt * sltaussamp
                    ADSLNMP(icol) = ADSLNMP(icol) + wgt * sltaussagmp
                    FORLDMP(icol) = FORLDMP(icol) + wgt * sltaussamp
                    FORLNMP(icol) = FORLNMP(icol) + wgt * sltaussafmp
                  end if
                  sitaussamp = 0.; sitaussagmp = 0.; sitaussafmp = 0.
                  if (sitaump > 0.) then
                    select type(cloud_props_gpt_ice)
                    class is (ty_optical_props_2str)
                      sitaussamp = sum(cloud_props_gpt_ice%tau(isub,LCLDMH:LCLDLM-1,igpt) * &
                                       cloud_props_gpt_ice%ssa(isub,LCLDMH:LCLDLM-1,igpt))
                      sitaussagmp = sum(cloud_props_gpt_ice%tau(isub,LCLDMH:LCLDLM-1,igpt) * &
                                        cloud_props_gpt_ice%ssa(isub,LCLDMH:LCLDLM-1,igpt) * &
                                        cloud_props_gpt_ice%g  (isub,LCLDMH:LCLDLM-1,igpt))
                      sitaussafmp = sum(cloud_props_gpt_ice%tau(isub,LCLDMH:LCLDLM-1,igpt) * &
                                        cloud_props_gpt_ice%ssa(isub,LCLDMH:LCLDLM-1,igpt) * &
                                                        forwice(isub,LCLDMH:LCLDLM-1,igpt))
                    end select
                    CDSIDMP(icol) = CDSIDMP(icol) + wgt
                    CDSINMP(icol) = CDSINMP(icol) + wgt * sitaump
                    SDSIDMP(icol) = SDSIDMP(icol) + wgt * sitaump
                    SDSINMP(icol) = SDSINMP(icol) + wgt * sitaussamp
                    ADSIDMP(icol) = ADSIDMP(icol) + wgt * sitaussamp
                    ADSINMP(icol) = ADSINMP(icol) + wgt * sitaussagmp
                    FORIDMP(icol) = FORIDMP(icol) + wgt * sitaussamp
                    FORINMP(icol) = FORINMP(icol) + wgt * sitaussafmp
                  end if

                  ! high pressure layer
                  sltauhp = sum(cloud_props_gpt_liq%tau(isub,1:LCLDMH-1,igpt))
                  sitauhp = sum(cloud_props_gpt_ice%tau(isub,1:LCLDMH-1,igpt))
                  stauhp = sltauhp + sitauhp
                  if (stauhp > 0.) then
                    CDSNHP(icol) = CDSNHP(icol) + wgt * stauhp
                    CDSDHP(icol) = CDSDHP(icol) + wgt
                  end if
                  sltaussahp = 0.; sltaussaghp = 0.; sltaussafhp = 0.
                  if (sltauhp > 0.) then
                    select type(cloud_props_gpt_liq)
                    class is (ty_optical_props_2str)
                      sltaussahp = sum(cloud_props_gpt_liq%tau(isub,1:LCLDMH-1,igpt) * &
                                       cloud_props_gpt_liq%ssa(isub,1:LCLDMH-1,igpt))
                      sltaussaghp = sum(cloud_props_gpt_liq%tau(isub,1:LCLDMH-1,igpt) * &
                                        cloud_props_gpt_liq%ssa(isub,1:LCLDMH-1,igpt) * &
                                        cloud_props_gpt_liq%g  (isub,1:LCLDMH-1,igpt))
                      sltaussafhp = sum(cloud_props_gpt_liq%tau(isub,1:LCLDMH-1,igpt) * &
                                        cloud_props_gpt_liq%ssa(isub,1:LCLDMH-1,igpt) * &
                                                        forwliq(isub,1:LCLDMH-1,igpt))
                    end select
                    CDSLDHP(icol) = CDSLDHP(icol) + wgt
                    CDSLNHP(icol) = CDSLNHP(icol) + wgt * sltauhp
                    SDSLDHP(icol) = SDSLDHP(icol) + wgt * sltauhp
                    SDSLNHP(icol) = SDSLNHP(icol) + wgt * sltaussahp
                    ADSLDHP(icol) = ADSLDHP(icol) + wgt * sltaussahp
                    ADSLNHP(icol) = ADSLNHP(icol) + wgt * sltaussaghp
                    FORLDHP(icol) = FORLDHP(icol) + wgt * sltaussahp
                    FORLNHP(icol) = FORLNHP(icol) + wgt * sltaussafhp
                  end if
                  sitaussahp = 0.; sitaussaghp = 0.; sitaussafhp = 0.
                  if (sitauhp > 0.) then
                    select type(cloud_props_gpt_ice)
                    class is (ty_optical_props_2str)
                      sitaussahp = sum(cloud_props_gpt_ice%tau(isub,1:LCLDMH-1,igpt) * &
                                       cloud_props_gpt_ice%ssa(isub,1:LCLDMH-1,igpt))
                      sitaussaghp = sum(cloud_props_gpt_ice%tau(isub,1:LCLDMH-1,igpt) * &
                                        cloud_props_gpt_ice%ssa(isub,1:LCLDMH-1,igpt) * &
                                        cloud_props_gpt_ice%g  (isub,1:LCLDMH-1,igpt))
                      sitaussafhp = sum(cloud_props_gpt_ice%tau(isub,1:LCLDMH-1,igpt) * &
                                        cloud_props_gpt_ice%ssa(isub,1:LCLDMH-1,igpt) * &
                                                        forwice(isub,1:LCLDMH-1,igpt))
                    end select
                    CDSIDHP(icol) = CDSIDHP(icol) + wgt
                    CDSINHP(icol) = CDSINHP(icol) + wgt * sitauhp
                    SDSIDHP(icol) = SDSIDHP(icol) + wgt * sitauhp
                    SDSINHP(icol) = SDSINHP(icol) + wgt * sitaussahp
                    ADSIDHP(icol) = ADSIDHP(icol) + wgt * sitaussahp
                    ADSINHP(icol) = ADSINHP(icol) + wgt * sitaussaghp
                    FORIDHP(icol) = FORIDHP(icol) + wgt * sitaussahp
                    FORINHP(icol) = FORINHP(icol) + wgt * sitaussafhp
                  end if

                  ! whole subcolumn
                  sltautp = sltaulp + sltaump + sltauhp
                  sitautp = sitaulp + sitaump + sitauhp
                  stautp = staulp + staump + stauhp
                  if (stautp > 0.) then
                    CDSNTP(icol) = CDSNTP(icol) + wgt * stautp
                    CDSDTP(icol) = CDSDTP(icol) + wgt
                  end if
                  sltaussatp = sltaussalp + sltaussamp + sltaussahp
                  sltaussagtp = sltaussaglp + sltaussagmp + sltaussaghp
                  sltaussaftp = sltaussaflp + sltaussafmp + sltaussafhp
                  if (sltautp > 0.) then
                    CDSLDTP(icol) = CDSLDTP(icol) + wgt
                    CDSLNTP(icol) = CDSLNTP(icol) + wgt * sltautp
                    SDSLDTP(icol) = SDSLDTP(icol) + wgt * sltautp
                    SDSLNTP(icol) = SDSLNTP(icol) + wgt * sltaussatp
                    ADSLDTP(icol) = ADSLDTP(icol) + wgt * sltaussatp
                    ADSLNTP(icol) = ADSLNTP(icol) + wgt * sltaussagtp
                    FORLDTP(icol) = FORLDTP(icol) + wgt * sltaussatp
                    FORLNTP(icol) = FORLNTP(icol) + wgt * sltaussaftp
                  end if
                  sitaussatp = sitaussalp + sitaussamp + sitaussahp
                  sitaussagtp = sitaussaglp + sitaussagmp + sitaussaghp
                  sitaussaftp = sitaussaflp + sitaussafmp + sitaussafhp
                  if (sitautp > 0.) then
                    CDSIDTP(icol) = CDSIDTP(icol) + wgt
                    CDSINTP(icol) = CDSINTP(icol) + wgt * sitautp
                    SDSIDTP(icol) = SDSIDTP(icol) + wgt * sitautp
                    SDSINTP(icol) = SDSINTP(icol) + wgt * sitaussatp
                    ADSIDTP(icol) = ADSIDTP(icol) + wgt * sitaussatp
                    ADSINTP(icol) = ADSINTP(icol) + wgt * sitaussagtp
                    FORIDTP(icol) = FORIDTP(icol) + wgt * sitaussatp
                    FORINTP(icol) = FORINTP(icol) + wgt * sitaussaftp
                  end if

                end do ! igpt
              end do ! ib

            end if  ! potentially cloudy column
          end do  ! isub
        end if  ! include_aerosols
        call MAPL_TimerOff(MAPL,"--RRTMGP_SPRLYR_DIAGS",__RC__)
#endif

        call MAPL_TimerOn(MAPL,"--RRTMGP_RT",__RC__)

        ! scale to our tsi
        ! (both toa_flux and tsi are NORMAL to solar beam, [W/m2])
        toa_flux = toa_flux * spread(tsi(colS:colE)/sum(toa_flux,dim=2), 2, ngpt)

        ! add in aerosol optical properties if requested and available
        if (need_aer_optical_props) then
          TEST_(aer_props%increment(optical_props))
        end if

        ! clear-sky radiative transfer
        fluxes_clrsky%flux_up  => flux_up_clrsky(colS:colE,:)
        fluxes_clrsky%flux_net => flux_net_clrsky(colS:colE,:)
        error_msg = rte_sw( &
          optical_props, top_at_1, mu0(colS:colE), toa_flux, &
          sfc_alb_dir(:,colS:colE), sfc_alb_dif(:,colS:colE), &
          fluxes_clrsky)
        TEST_(error_msg)

        ! add in cloud optical properties
        ! add ice first since its optical depths are usually smaller
        TEST_(cloud_props_gpt_ice%increment(optical_props))
        TEST_(cloud_props_gpt_liq%increment(optical_props))

        ! all-sky radiative transfer
        fluxes_allsky%flux_up         => flux_up_allsky(colS:colE,:)
        fluxes_allsky%flux_net        => flux_net_allsky(colS:colE,:)
        fluxes_allsky%bnd_flux_dn     => bnd_flux_dn_allsky(colS:colE,:,:)
        fluxes_allsky%bnd_flux_dn_dir => bnd_flux_dir_allsky(colS:colE,:,:)
        fluxes_allsky%bnd_flux_net    => bnd_flux_net_allsky(colS:colE,:,:)
        error_msg = rte_sw( &
          optical_props, top_at_1, mu0(colS:colE), toa_flux, &
          sfc_alb_dir(:,colS:colE), sfc_alb_dif(:,colS:colE), &
          fluxes_allsky)
        TEST_(error_msg)

        call MAPL_TimerOff(MAPL,"--RRTMGP_RT",__RC__)

      end do ! loop over blocks

      call MAPL_TimerOn(MAPL,"--RRTMGP_POST",__RC__)

      ! normalize by incoming solar radiation
      allocate(flux_dn_top(ncol), __STAT__)
      flux_dn_top(:) = real(max(SLR1D, 1e-7), kind=wp)
      do k = 1, LM+1
        flux_up_clrsky (:,k) = flux_up_clrsky (:,k) / flux_dn_top(:)
        flux_net_clrsky(:,k) = flux_net_clrsky(:,k) / flux_dn_top(:)
      end do
      do k = 1, LM+1
        flux_up_allsky (:,k) = flux_up_allsky (:,k) / flux_dn_top(:)
        flux_net_allsky(:,k) = flux_net_allsky(:,k) / flux_dn_top(:)
      end do
      do ib = 1, nbnd
        do k = 1, LM+1
          bnd_flux_dn_allsky (:,k,ib) = bnd_flux_dn_allsky (:,k,ib) / flux_dn_top(:)
          bnd_flux_net_allsky(:,k,ib) = bnd_flux_net_allsky(:,k,ib) / flux_dn_top(:)
          bnd_flux_dir_allsky(:,k,ib) = bnd_flux_dir_allsky(:,k,ib) / flux_dn_top(:)
        end do
      end do
      deallocate(flux_dn_top, __STAT__)

      ! load output arrays
      ! clear-sky fluxes
      FSCU = real(flux_up_clrsky)
      FSC  = real(flux_net_clrsky)
      ! all-sky fluxes
      FSWU = real(flux_up_allsky)
      FSW  = real(flux_net_allsky)
      ! surface net flux per band
      do ib = 1, nbnd
        FSWBAND(:,ib) = real(bnd_flux_net_allsky(:,LM+1,ib))
      end do
      ! surface downwelling direct and diffuse fluxes in bands
      if (SOLAR_TO_OBIO .and. include_aerosols) then
         do ib = 1, nbnd
            DRBAND(:,ib) = real(bnd_flux_dir_allsky(:,LM+1,ib))
            DFBAND(:,ib) = real(bnd_flux_dn_allsky (:,LM+1,ib) - bnd_flux_dir_allsky(:,LM+1,ib))
         end do
      endif

      ! surface direct and diffuse downward in super-bands
      ! for *diffuse* downward must subtract direct (downward) from total downward
      ! pmn: may later do this using a flux class extension??

      ! NIR bands (1-9: 820-12850 cm-1, 0.778-12.195 microns)
      NIRR = 0.; NIRF = 0.
      do ib=1,9
        NIRR = NIRR + real(bnd_flux_dir_allsky(:,LM+1,ib))
        NIRF = NIRF + real(bnd_flux_dn_allsky (:,LM+1,ib) - bnd_flux_dir_allsky(:,LM+1,ib))
      end do
      ! PAR bands (11-12: 16000-29000 cm-1, 0.345-0.625 micron)
      PARR = 0.; PARF = 0.
      do ib=11,12
        PARR = PARR + real(bnd_flux_dir_allsky(:,LM+1,ib))
        PARF = PARF + real(bnd_flux_dn_allsky (:,LM+1,ib) - bnd_flux_dir_allsky(:,LM+1,ib))
      enddo
      ! UVR bands (13-14: 29000-50000 cm-1, 0.200-0.345 micron)
      UVRR = 0.; UVRF = 0.
      do ib=13,14
        UVRR = UVRR + real(bnd_flux_dir_allsky(:,LM+1,ib))
        UVRF = UVRF + real(bnd_flux_dn_allsky (:,LM+1,ib) - bnd_flux_dir_allsky(:,LM+1,ib))
      enddo
      ! Transition band (10, 12850-16000 cm-1, 0.625-0.778 micron)
      ! split half-and-half to PAR and NIR
      NIRR = NIRR + 0.5 * real(bnd_flux_dir_allsky(:,LM+1,10))
      PARR = PARR + 0.5 * real(bnd_flux_dir_allsky(:,LM+1,10))
      NIRF = NIRF + 0.5 * real(bnd_flux_dn_allsky (:,LM+1,10) - bnd_flux_dir_allsky(:,LM+1,10))
      PARF = PARF + 0.5 * real(bnd_flux_dn_allsky (:,LM+1,10) - bnd_flux_dir_allsky(:,LM+1,10))

      ! clean up
      deallocate(band_lims_gpt,__STAT__)
      deallocate(tsi,mu0,sfc_alb_dir,sfc_alb_dif,toa_flux,__STAT__)
      deallocate(dummy_wp,p_lay,t_lay,p_lev,dp_wp,dzmid,__STAT__)
      deallocate(flux_up_clrsky,flux_net_clrsky,__STAT__)
      deallocate(flux_up_allsky,flux_net_allsky,__STAT__)
      deallocate(bnd_flux_dn_allsky,bnd_flux_net_allsky,bnd_flux_dir_allsky,__STAT__)
      deallocate(seeds,urand,cld_mask,__STAT__)
      deallocate(forwliq,forwice,__STAT__)
      if (gen_mro) then
        deallocate(adl,alpha,urand_aux,__STAT__)
        if (cond_inhomo) then
          deallocate(rdl,rcorr,urand_cond,urand_cond_aux,zcw,__STAT__)
        endif
      end if
      if (include_aerosols) deallocate(ClearCounts,__STAT__)
      call cloud_optics%finalize()
      call cloud_props_gpt_liq%finalize()
      call cloud_props_gpt_ice%finalize()
      call cloud_props_bnd_liq%finalize()
      call cloud_props_bnd_ice%finalize()
      if (need_aer_optical_props) call aer_props%finalize()
      call optical_props%finalize()

      call MAPL_TimerOff(MAPL,"--RRTMGP_POST",__RC__)

      call MAPL_TimerOff(MAPL,"-RRTMGP",__RC__)

#undef TEST_

   else if (USE_RRTMG) then

      ! regular RRTMG
      call MAPL_TimerOn(MAPL,"-RRTMG")

      NCOL = size(Q,1)

      ! reversed (flipped) vertical dimension arrays and other RRTMG arrays
      ! -------------------------------------------------------------------

      ! interface (between layer) variables
      allocate(TLEV  (size(Q,1),size(Q,2)+1),__STAT__)
      allocate(TLEV_R(size(Q,1),size(Q,2)+1),__STAT__)
      allocate(PLE_R (size(Q,1),size(Q,2)+1),__STAT__)
      ! cloud physical properties
      allocate(FCLD_R(size(Q,1),size(Q,2)),__STAT__)
      allocate(CLIQWP(size(Q,1),size(Q,2)),__STAT__)
      allocate(CICEWP(size(Q,1),size(Q,2)),__STAT__)
      allocate(RELIQ (size(Q,1),size(Q,2)),__STAT__)
      allocate(REICE (size(Q,1),size(Q,2)),__STAT__)
      ! aerosol optical properties
      allocate(TAUAER(size(Q,1),size(Q,2),NB_RRTMG),__STAT__)
      allocate(SSAAER(size(Q,1),size(Q,2),NB_RRTMG),__STAT__)
      allocate(ASMAER(size(Q,1),size(Q,2),NB_RRTMG),__STAT__)
      ! layer variables
      allocate(DPR   (size(Q,1),size(Q,2)),__STAT__)
      allocate(PL_R  (size(Q,1),size(Q,2)),__STAT__)
      allocate(ZL_R  (size(Q,1),size(Q,2)),__STAT__)
      allocate(T_R   (size(Q,1),size(Q,2)),__STAT__)
      allocate(Q_R   (size(Q,1),size(Q,2)),__STAT__)
      allocate(O2_R  (size(Q,1),size(Q,2)),__STAT__)
      allocate(O3_R  (size(Q,1),size(Q,2)),__STAT__)
      allocate(CO2_R (size(Q,1),size(Q,2)),__STAT__)
      allocate(CH4_R (size(Q,1),size(Q,2)),__STAT__)
      ! super-layer cloud fractions
      allocate(CLEARCOUNTS (size(Q,1),4),__STAT__)
      ! output fluxes
      allocate(SWUFLX (size(Q,1),size(Q,2)+1),__STAT__)
      allocate(SWDFLX (size(Q,1),size(Q,2)+1),__STAT__)
      allocate(SWUFLXC(size(Q,1),size(Q,2)+1),__STAT__)
      allocate(SWDFLXC(size(Q,1),size(Q,2)+1),__STAT__)
      ! un-flipped outputs
      allocate(SWUFLXR (size(Q,1),size(Q,2)+1),__STAT__)
      allocate(SWDFLXR (size(Q,1),size(Q,2)+1),__STAT__)
      allocate(SWUFLXCR(size(Q,1),size(Q,2)+1),__STAT__)
      allocate(SWDFLXCR(size(Q,1),size(Q,2)+1),__STAT__)

      ! Set flags related to cloud properties (see RRTMG_SW)
      ! ----------------------------------------------------
      call MAPL_GetResource(MAPL,ICEFLGSW,'RRTMG_ICEFLG:',DEFAULT=3,__RC__)
      call MAPL_GetResource(MAPL,LIQFLGSW,'RRTMG_LIQFLG:',DEFAULT=1,__RC__)

      ! Normalize aerosol inputs
      ! ------------------------

      if (num_aero_vars > 0) then
         where (TAUA > 0. .and. SSAA > 0. )
            ASYA = ASYA/SSAA
            SSAA = SSAA/TAUA
         elsewhere
            TAUA = 0.
            SSAA = 0.
            ASYA = 0.
         end where
      end if

      ! Flip in vertical, Convert units, and interpolate T, etc.
      ! --------------------------------------------------------
      ! RRTMG convention is that vertical indices increase from bot -> top

      call MAPL_TimerOn (MAPL,"--RRTMG_FLIP")

      DPR(:,1:LM) = (PLE(:,2:LM+1)-PLE(:,1:LM))

      ! cloud water paths converted from g/g to g/m^2
      CICEWP(:,1:LM) = (1.02*100*DPR(:,LM:1:-1))*QQ3(:,LM:1:-1,1)
      CLIQWP(:,1:LM) = (1.02*100*DPR(:,LM:1:-1))*QQ3(:,LM:1:-1,2)

      ! cloud effective radii with limits imposed as assumed by RRTMG
      REICE (:,1:LM) = RR3(:,LM:1:-1,1)
      RELIQ (:,1:LM) = RR3(:,LM:1:-1,2)

      IF      (ICEFLGSW == 0) THEN
         WHERE (REICE < 10.)  REICE = 10.
         WHERE (REICE > 30.)  REICE = 30.
      ELSE IF (ICEFLGSW == 1) THEN
         WHERE (REICE < 13.)  REICE = 13.
         WHERE (REICE > 130.) REICE = 130.
      ELSE IF (ICEFLGSW == 2) THEN
         WHERE (REICE < 5.)   REICE = 5.
         WHERE (REICE > 131.) REICE = 131.
      ELSE IF (ICEFLGSW == 3) THEN
         WHERE (REICE < 5.)   REICE = 5.
         WHERE (REICE > 140.) REICE = 140.
      ELSE IF (ICEFLGSW == 4) THEN
         REICE(:,:) = REICE(:,:)*2.
         WHERE (REICE < 1.)   REICE = 1.
         WHERE (REICE > 200.) REICE = 200.
      END IF

      IF      (LIQFLGSW == 0) THEN
         WHERE (RELIQ < 10.)  RELIQ = 10.
         WHERE (RELIQ > 30.)  RELIQ = 30.
      ELSE IF (LIQFLGSW == 1) THEN
         WHERE (RELIQ < 2.5)  RELIQ = 2.5
         WHERE (RELIQ > 60.)  RELIQ = 60.
      END IF

      ! regular (non-flipped) interface temperatures
      TLEV(:,2:LM)=(T(:,1:LM-1)* DPR(:,2:LM) + T(:,2:LM) * DPR(:,1:LM-1)) &
            /(DPR(:,1:LM-1) + DPR(:,2:LM))
      TLEV(:,LM+1) = TS(:)
      TLEV(:,   1) = TLEV(:,2)

      ! flip in vertical ...

      PLE_R (:,1:LM+1) = PLE (:,LM+1:1:-1) / 100.  ! hPa
      TLEV_R(:,1:LM+1) = TLEV(:,LM+1:1:-1)

      PL_R  (:,1:LM  ) = PL  (:,LM:1:-1)   / 100.  ! hPa
      T_R   (:,1:LM  ) = T   (:,LM:1:-1)

      ! Specific humidity is converted to Volume Mixing Ratio
      Q_R   (:,1:LM  ) = Q  (:,LM:1:-1) / (1.-Q(:,LM:1:-1)) * (MAPL_AIRMW/MAPL_H2OMW)

      ! Ozone is converted Mass Mixing Ratio to Volume Mixing Ratio
      O3_R  (:,1:LM  ) = O3 (:,LM:1:-1) * (MAPL_AIRMW/MAPL_O3MW)

      ! chemistry and cloud fraction
      ! (cloud water paths and effective radii flipped already)
      CH4_R (:,1:LM  ) = CH4(:,LM:1:-1)
      CO2_R (:,1:LM  ) = CO2
      O2_R  (:,1:LM  ) = O2
      FCLD_R(:,1:LM  ) = CL (:,LM:1:-1)

! Clean up negatives
      WHERE (Q_R < 0.) Q_R = 0.
      WHERE (O3_R < 0.) O3_R = 0.
      WHERE (CH4_R < 0.) CH4_R = 0.
      WHERE (CO2_R < 0.) CO2_R = 0.
      WHERE (O2_R < 0.) O2_R = 0.
      WHERE (FCLD_R < 0.) FCLD_R = 0.

      ! Adjustment for Earth/Sun distance, from MAPL_SunGetInsolation
      ADJES = DIST

      ! Layer mid-point heights relative to zero at index 1
      ZL_R(:,1) = 0.
      do k=2,LM
         ! dz = RT/g x dp/p
         ! Note: This is correct even though its different from LW.
         ! Its because SW uses LE[V]_R 1:LM+1 while LW uses 0:LM.
         ZL_R(:,k) = ZL_R(:,k-1) + MAPL_RGAS*TLEV_R(:,k)/MAPL_GRAV*(PL_R(:,k-1)-PL_R(:,k))/PLE_R(:,k)
      enddo

      ! aerosols
      TAUAER(:,1:LM,:) = TAUA(:,LM:1:-1,:)
      SSAAER(:,1:LM,:) = SSAA(:,LM:1:-1,:)
      ASMAER(:,1:LM,:) = ASYA(:,LM:1:-1,:)

      call MAPL_TimerOff(MAPL,"--RRTMG_FLIP")
      call MAPL_TimerOn (MAPL,"--RRTMG_INIT")

      ! initialize RRTMG SW
      call RRTMG_SW_INI

      call MAPL_TimerOff(MAPL,"--RRTMG_INIT")
      call MAPL_TimerOn (MAPL,"--RRTMG_RUN")

      ! partition size for columns (profiles) used to improve efficiency
      call MAPL_GetResource( MAPL, RPART, 'RRTMGSW_PARTITION_SIZE:',  DEFAULT=0, __RC__)

      ! various RRTMG configuration options ...

      IAER = 10    ! Per AER:
                   !  0: Turns off aerosols
                   ! 10: Enables aerosols

      NORMFLX = 1  ! 0: Do not normalize fluxes
                   ! 1: Normalize fluxes

      DYOFYR = DOY ! Day of year

      call MAPL_GetResource( MAPL, ISOLVAR ,'ISOLVAR:', DEFAULT=0, __RC__)

                   ! ISOLVAR:
                   ! Flag for solar variability method
                   !   -1 = (when scon .eq. 0.0): No solar variability
                   !        and no solar cycle (Kurucz solar irradiance
                   !        of 1368.22 Wm-2 only);
                   !        (when scon .ne. 0.0): Kurucz solar irradiance
                   !        scaled to scon and solar variability defined
                   !        (optional) by setting non-zero scale factors
                   !        for each band in bndsolvar
                   !    0 = (when SCON .eq. 0.0): No solar variability
                   !        and no solar cycle (NRLSSI2 solar constant of
                   !        1360.85 Wm-2 for the 100-50000 cm-1 spectral
                   !        range only), with facular and sunspot effects
                   !        fixed to the mean of Solar Cycles 13-24;
                   !        (when SCON .ne. 0.0): No solar variability
                   !        and no solar cycle (NRLSSI2 solar constant of
                   !        1360.85 Wm-2 for the 100-50000 cm-1 spectral
                   !        range only), is scaled to SCON
                   !    1 = Solar variability (using NRLSSI2  solar
                   !        model) with solar cycle contribution
                   !        determined by fraction of solar cycle
                   !        with facular and sunspot variations
                   !        fixed to their mean variations over the
                   !        average of Solar Cycles 13-24;
                   !        two amplitude scale factors allow
                   !        facular and sunspot adjustments from
                   !        mean solar cycle as defined by indsolvar
                   !    2 = Solar variability (using NRLSSI2 solar
                   !        model) over solar cycle determined by
                   !        direct specification of Mg (facular)
                   !        and SB (sunspot) indices provided
                   !        in indsolvar (scon = 0.0 only)
                   !    3 = (when scon .eq. 0.0): No solar variability
                   !        and no solar cycle (NRLSSI2 solar irradiance
                   !        of 1360.85 Wm-2 only);
                   !        (when scon .ne. 0.0): NRLSSI2 solar irradiance
                   !        scaled to scon and solar variability defined
                   !        (optional) by setting non-zero scale factors
                   !        for each band in bndsolvar

      if (ISOLVAR == 1) then
         if (MAPL_AM_I_ROOT()) then
            write (*,*) "ERROR in RRTMG_SW"
            write (*,*) "ISOLVAR==1 is currently unsupported as we have no"
            write (*,*) "way of correctly setting solcycfrac."
         end if
         _FAIL('RRTMG SW: ISOLVAR==1 currently unsupported')
      end if

      ! INDSOLVAR =  Facular and sunspot amplitude
      !              scale factors (isolvar=1), or
      !              Mg and SB indices (isolvar=2)
      !                 Dimensions: (2)

      if (ISOLVAR == 2) then

         ! Solar indices from our file
         INDSOLVAR(1) = MG
         INDSOLVAR(2) = SB

      else

         call MAPL_GetResource( MAPL, INDSOLVAR(1) ,'INDSOLVAR_1:', DEFAULT=1.0, __RC__)
         call MAPL_GetResource( MAPL, INDSOLVAR(2) ,'INDSOLVAR_2:', DEFAULT=1.0, __RC__)

      end if


      BNDSOLVAR = 1.0 ! Solar variability scale factors
                      ! for each shortwave band
                      !    Dimensions: (nbndsw=14)

      ! SOLCYCFRAC: Fraction of averaged 11-year solar cycle (0-1)
      !                at current time (isolvar=1)
      !                0.0 represents the first day of year 1
      !                1.0 represents the last day of year 11

      ! MAT: Note while we don't currently use SOLCYCFRAC, we set it to something
      !      to avoid an optional variable on GPUs

      call MAPL_GetResource( MAPL, SOLCYCFRAC ,'SOLCYCFRAC:', DEFAULT=1.0, __RC__)

      ! call RRTMG SW
      ! -------------

      call RRTMG_SW (MAPL, &
         RPART, NCOL, LM, &
         SC, ADJES, ZT, ISOLVAR, &
         PL_R, PLE_R, T_R, &
         Q_R, O3_R, CO2_R, CH4_R, O2_R, &
         ICEFLGSW, LIQFLGSW, &
         FCLD_R, CICEWP, CLIQWP, REICE, RELIQ, &
         DYOFYR, ZL_R, ALAT, &
         IAER, TAUAER, SSAAER, ASMAER, &
         ALBVR, ALBVF, ALBNR, ALBNF, &
         LM-LCLDLM+1, LM-LCLDMH+1, NORMFLX, &
         CLEARCOUNTS, SWUFLX, SWDFLX, SWUFLXC, SWDFLXC, &
         NIRR, NIRF, PARR, PARF, UVRR, UVRF, FSWBAND, &

         COTDTP, COTDHP, COTDMP, COTDLP, &
         COTNTP, COTNHP, COTNMP, COTNLP, &

#ifdef SOLAR_RADVAL
         CDSDTP, CDSDHP, CDSDMP, CDSDLP, &
         CDSNTP, CDSNHP, CDSNMP, CDSNLP, &

         COTLDTP, COTLDHP, COTLDMP, COTLDLP, &
         COTLNTP, COTLNHP, COTLNMP, COTLNLP, &
         CDSLDTP, CDSLDHP, CDSLDMP, CDSLDLP, &
         CDSLNTP, CDSLNHP, CDSLNMP, CDSLNLP, &
         COTIDTP, COTIDHP, COTIDMP, COTIDLP, &
         COTINTP, COTINHP, COTINMP, COTINLP, &
         CDSIDTP, CDSIDHP, CDSIDMP, CDSIDLP, &
         CDSINTP, CDSINHP, CDSINMP, CDSINLP, &

         SSALDTP, SSALDHP, SSALDMP, SSALDLP, &
         SSALNTP, SSALNHP, SSALNMP, SSALNLP, &
         SDSLDTP, SDSLDHP, SDSLDMP, SDSLDLP, &
         SDSLNTP, SDSLNHP, SDSLNMP, SDSLNLP, &
         SSAIDTP, SSAIDHP, SSAIDMP, SSAIDLP, &
         SSAINTP, SSAINHP, SSAINMP, SSAINLP, &
         SDSIDTP, SDSIDHP, SDSIDMP, SDSIDLP, &
         SDSINTP, SDSINHP, SDSINMP, SDSINLP, &

         ASMLDTP, ASMLDHP, ASMLDMP, ASMLDLP, &
         ASMLNTP, ASMLNHP, ASMLNMP, ASMLNLP, &
         ADSLDTP, ADSLDHP, ADSLDMP, ADSLDLP, &
         ADSLNTP, ADSLNHP, ADSLNMP, ADSLNLP, &
         ASMIDTP, ASMIDHP, ASMIDMP, ASMIDLP, &
         ASMINTP, ASMINHP, ASMINMP, ASMINLP, &
         ADSIDTP, ADSIDHP, ADSIDMP, ADSIDLP, &
         ADSINTP, ADSINHP, ADSINMP, ADSINLP, &

         FORLDTP, FORLDHP, FORLDMP, FORLDLP, &
         FORLNTP, FORLNHP, FORLNMP, FORLNLP, &
         FORIDTP, FORIDHP, FORIDMP, FORIDLP, &
         FORINTP, FORINHP, FORINMP, FORINLP, &
#endif

         SOLAR_TO_OBIO .and. include_aerosols, DRBAND, DFBAND, &
         BNDSOLVAR, INDSOLVAR, SOLCYCFRAC, &
         __RC__)

      call MAPL_TimerOff(MAPL,"--RRTMG_RUN")
      call MAPL_TimerOn (MAPL,"--RRTMG_FLIP")

      ! unflip the outputs in the vertical
      ! ----------------------------------

      SWUFLXR (:,1:LM+1) = SWUFLX (:,LM+1:1:-1)
      SWDFLXR (:,1:LM+1) = SWDFLX (:,LM+1:1:-1)
      SWUFLXCR(:,1:LM+1) = SWUFLXC(:,LM+1:1:-1)
      SWDFLXCR(:,1:LM+1) = SWDFLXC(:,LM+1:1:-1)

      call MAPL_TimerOff(MAPL,"--RRTMG_FLIP")

      ! required outputs
      ! ----------------

      ! convert super-layer clearCounts to cloud fractions
      if (include_aerosols) then
        CLDTS(:) = 1. - CLEARCOUNTS(:,1)/float(NGPTSW)
        CLDHS(:) = 1. - CLEARCOUNTS(:,2)/float(NGPTSW)
        CLDMS(:) = 1. - CLEARCOUNTS(:,3)/float(NGPTSW)
        CLDLS(:) = 1. - CLEARCOUNTS(:,4)/float(NGPTSW)
      end if

      ! undef versions of cloud optical thicknesses
      ! We cannot use a merge here because some compilers (e.g. Intel)
      ! will will evaluate the first argument first and if both are
      ! zero, it will return NaNs.
      where (COTNTP > 0. .and. COTDTP > 0.)
         COTTP = COTNTP/COTDTP
      elsewhere
         COTTP = MAPL_UNDEF
      end where

      where (COTNHP > 0. .and. COTDHP > 0.)
         COTHP = COTNHP/COTDHP
      elsewhere
         COTHP = MAPL_UNDEF
      end where

      where (COTNMP > 0. .and. COTDMP > 0.)
         COTMP = COTNMP/COTDMP
      elsewhere
         COTMP = MAPL_UNDEF
      end where

      where (COTNLP > 0. .and. COTDLP > 0.)
         COTLP = COTNLP/COTDLP
      elsewhere
         COTLP = MAPL_UNDEF
      end where

#ifdef SOLAR_RADVAL
      ! zero versions of cloud optical thicknesses
      ! We can use merge() here because we cannot divide by zero
      TAUTP = merge(COTTP, 0., COTNTP > 0. .and. COTDTP > 0.)
      TAUHP = merge(COTHP, 0., COTNHP > 0. .and. COTDHP > 0.)
      TAUMP = merge(COTMP, 0., COTNMP > 0. .and. COTDMP > 0.)
      TAULP = merge(COTLP, 0., COTNLP > 0. .and. COTDLP > 0.)
#endif

      ! fluxes
      FSW  = SWDFLXR  - SWUFLXR
      FSC  = SWDFLXCR - SWUFLXCR
      FSWU = SWUFLXR
      FSCU = SWUFLXCR

      ! Deallocate the working inputs
      !------------------------------
      deallocate(TLEV  ,__STAT__)
      deallocate(TLEV_R,__STAT__)
      deallocate(PLE_R ,__STAT__)
      deallocate(FCLD_R,__STAT__)
      deallocate(CLIQWP,__STAT__)
      deallocate(CICEWP,__STAT__)
      deallocate(RELIQ ,__STAT__)
      deallocate(REICE ,__STAT__)

      deallocate(TAUAER,__STAT__)
      deallocate(SSAAER,__STAT__)
      deallocate(ASMAER,__STAT__)
      deallocate(DPR   ,__STAT__)
      deallocate(PL_R  ,__STAT__)
      deallocate(ZL_R  ,__STAT__)
      deallocate(T_R   ,__STAT__)
      deallocate(Q_R   ,__STAT__)
      deallocate(O2_R  ,__STAT__)
      deallocate(O3_R  ,__STAT__)
      deallocate(CO2_R ,__STAT__)
      deallocate(CH4_R ,__STAT__)

      deallocate(CLEARCOUNTS ,__STAT__)

      deallocate(SWUFLX ,__STAT__)
      deallocate(SWDFLX ,__STAT__)
      deallocate(SWUFLXC,__STAT__)
      deallocate(SWDFLXC,__STAT__)

      deallocate(SWUFLXR ,__STAT__)
      deallocate(SWDFLXR ,__STAT__)
      deallocate(SWUFLXCR,__STAT__)
      deallocate(SWDFLXCR,__STAT__)

      call MAPL_TimerOff(MAPL,"-RRTMG")

   else

      _FAIL('unknown SW radiation scheme!')

   end if SCHEME

      ! Deallocate the working inputs
      !------------------------------

      deallocate (PL, RH, PLhPa)
      deallocate (QQ3, RR3)
      deallocate (O3)
      deallocate (TAUA, SSAA, ASYA)

      ! Complete load balancing by retrieving work done remotely
      !---------------------------------------------------------

      call MAPL_TimerOn(MAPL,"-BALANCE")

      call MAPL_TimerOn(MAPL,"--RETRIEVE")
      if (LoadBalance) then
         if (size(BufOut)   > 0) call MAPL_BalanceWork(BufOut,  NumMax,Direction=MAPL_Retrieve,Handle=SolarBalanceHandle,__RC__)
         if (size(BufInOut) > 0) call MAPL_BalanceWork(BufInOut,NumMax,Direction=MAPL_Retrieve,Handle=SolarBalanceHandle,__RC__)
      end if
      call MAPL_TimerOff(MAPL,"--RETRIEVE")

      ! Unpack the results. Fills "masked" (night) locations with default value from internal state
      !--------------------------------------------------------------------------------------------
      ! resulting internals are then contiguous versions
      ! Note: InOut variables do not fill unmasked locations with a default,
      ! since the unmasked locations may contain potentially useful aged data.

      i1InOut = 1; i1Out = 1
      INT_VARS_3: do k=1,NumInt
         if (SlicesInt(k) == 0) cycle

         if (IntInOut(k)) then
            pi1 => i1InOut
         else
            pi1 => i1Out
            call MAPL_VarSpecGet(InternalSpec(k),DEFAULT=def,__RC__)
         endif

         if (ugDim(k) > 0) then
            select case(rgDim(k))
               case(MAPL_DIMSHORZVERT)
                  call ESMFL_StateGetPointerToData(INTERNAL,ptr4,NamesInt(k),__RC__)
                  if (IntInOut(k)) then
                     do j=1,ugDim(k)
                        call UnPackIt(BufInOut(pi1+(j-1)*size(ptr4,3)*NumMax),ptr4(:,:,:,j), &
                           daytime,NumMax,HorzDims,size(ptr4,3))
                     end do
                  else
                     do j=1,ugDim(k)
                        call UnPackIt(BufOut  (pi1+(j-1)*size(ptr4,3)*NumMax),ptr4(:,:,:,j), &
                           daytime,NumMax,HorzDims,size(ptr4,3),def)
                     end do
                  endif
               case(MAPL_DIMSHORZONLY)
                  call ESMFL_StateGetPointerToData(INTERNAL,ptr3,NamesInt(k),__RC__)
                  if (IntInOut(k)) then
                     call UnPackIt(BufInOut(pi1),ptr3,daytime,NumMax,HorzDims,ugDim(k))
                  else
                     call UnPackIt(BufOut  (pi1),ptr3,daytime,NumMax,HorzDims,ugDim(k),def)
                  endif
            end select
         else
            select case(rgDim(k))
               case(MAPL_DIMSHORZVERT)
                  call ESMFL_StateGetPointerToData(INTERNAL,ptr3,NamesInt(k),__RC__)
                  if (IntInOut(k)) then
                     call UnPackIt(BufInOut(pi1),ptr3,daytime,NumMax,HorzDims,size(ptr3,3))
                  else
                     call UnPackIt(BufOut  (pi1),ptr3,daytime,NumMax,HorzDims,size(ptr3,3),def)
                  endif
               case(MAPL_DIMSHORZONLY)
                  call ESMFL_StateGetPointerToData(INTERNAL,ptr2,NamesInt(k),__RC__)
                  if (IntInOut(k)) then
                     call UnPackIt(BufInOut(pi1),ptr2,daytime,NumMax,HorzDims,1)
                  else
                     call UnPackIt(BufOut  (pi1),ptr2,daytime,NumMax,HorzDims,1,def)
                  endif
            end select
         end if
         pi1 = pi1 + NumMax*SlicesInt(k)

      enddo INT_VARS_3

      ! clean up
      deallocate(SlicesInp,NamesInp,__STAT__)
      deallocate(SlicesInt,NamesInt,__STAT__)
      deallocate(IntInOut,rgDim,ugDim,__STAT__)
      deallocate(BufInp,BufInOut,BufOut,__STAT__)
      call MAPL_TimerOn(MAPL,"--DESTROY")
      if (LoadBalance) call MAPL_BalanceDestroy(Handle=SolarBalanceHandle, __RC__)
      call MAPL_TimerOff(MAPL,"--DESTROY")

      call MAPL_TimerOff(MAPL,"-BALANCE")

      RETURN_(ESMF_SUCCESS)
    end subroutine SORADCORE


    subroutine SHRTWAVE(PLhPa,TA,WA,OA,CO2,COSZ   , &
                          CWC,REFF,FCLD,ICT,ICB  , &
                          RGBUV,RGFUV,RGBIR,RGFIR, &
                          TAUA,SSAA,ASYA, &

                          FLX,    FLC     , &
                          FDIRIR, FDIFIR  , &
                          FDIRPAR,FDIFPAR , &
                          FDIRUV, FDIFUV  , &
                          FLXU,   FLCU    , &
                          FLXBAND,          &
                          do_drfband,       &
                          drband, dfband,   &

                          RC                )


!   Inlineable cover for the f77 version of SORAD.
!   This cover works on a 1D run of soundings.

      real, dimension(:,:    ), intent(IN ) :: PLhPa,TA, WA, OA, FCLD
      real, dimension(:,:,:  ), intent(IN ) :: CWC, REFF
      real,                     intent(IN ) :: CO2
      integer,                  intent(IN ) :: ICT, ICB
      real, dimension(:      ), intent(IN ) :: COSZ, RGBUV, RGFUV, RGBIR, RGFIR
      real, dimension(:,:,:  ), intent(IN ) :: TAUA, SSAA, ASYA
      logical,                  intent(in ) :: do_drfband  ! Compute drband, dfband?

      real, dimension(:,:    ), intent(OUT) :: FLX,    FLC
      real, dimension(:,:    ), intent(OUT) :: FLXU,   FLCU
      real, dimension(:      ), intent(OUT) :: FDIRIR, FDIFIR
      real, dimension(:      ), intent(OUT) :: FDIRPAR,FDIFPAR
      real, dimension(:      ), intent(OUT) :: FDIRUV, FDIFUV
      real, dimension(:,:    ), intent(OUT) :: FLXBAND

      ! Surface downwelling direct and diffuse (W/m2) in each solar band:
      ! Only filled if (do_drfband), otherwise not touched and can be null pointers;
      ! if (do_drfband), must point to an (#cols,#bands) space.
      real, intent(inout), dimension (:,:), pointer :: drband, dfband

      integer, optional,        intent(OUT) :: RC

! Locals
!-------

      integer :: IRUN, LN
      INTEGER :: STATUS

! Begin
!------

      call MAPL_TimerOn(MAPL,"-MISC")

      IRUN     = SIZE(TA,1)
      LN       = SIZE(TA,2)

      call MAPL_TimerOff(MAPL,"-MISC")

      call MAPL_TimerOn(MAPL,"-SORAD")

      call MAPL_TimerOn(MAPL,"--SORAD_RUN",__RC__)
      call SORAD (IRUN,LN,NB_CHOU,COSZ,PLhPa,TA,WA,OA,CO2,      &
           CWC,FCLD,ICT,ICB,REFF,HK_UV_TEMP,HK_IR_TEMP,         &
           TAUA,SSAA,ASYA,                                      &
           RGBUV,RGFUV, RGBIR, RGFIR,                           &
           FLX,FLC,FDIRUV,FDIFUV,FDIRPAR,FDIFPAR,FDIRIR,FDIFIR, &
           FLXU,FLCU,                                           &
           FLXBAND,                                             &
           do_drfband, drband, dfband)
      call MAPL_TimerOff(MAPL,"--SORAD_RUN",__RC__)

      call MAPL_TimerOff(MAPL,"-SORAD")

      RETURN_(ESMF_SUCCESS)

    end subroutine SHRTWAVE

!==========================================================

    subroutine UPDATE_EXPORT(IM,JM,LM, RC)

      use parrrsw, only: nbndsw, jpb1, jpb2
      use rrsw_wvn, only: wavenum1, wavenum2
      use mo_gas_concentrations, only: ty_gas_concs
      use mo_load_coefficients,  only: load_and_init

      integer,           intent(IN ) :: IM, JM, LM
      integer, optional, intent(OUT) :: RC

!  Locals

      character(len=ESMF_MAXSTR)      :: IAm
      integer                         :: STATUS

      real,    dimension(IM,JM)       :: ZTH, SLR, ALB, CLD, SLN, ZTHN

      real, pointer, dimension(:,:  ) :: ALBEXP, ALBIMP

      real, pointer, dimension(:,:,:) ::     FSW,     FSC
      real, pointer, dimension(:,:,:) ::    FSWN,    FSCN
      real, pointer, dimension(:,:,:) ::   FSWNA,   FSCNA
      real, pointer, dimension(:,:,:) ::  FSWNAN,  FSCNAN

      real, pointer, dimension(:,:,:) ::    FSWU,    FSCU
      real, pointer, dimension(:,:,:) ::   FSWUN,   FSCUN
      real, pointer, dimension(:,:,:) ::  FSWUNA,  FSCUNA
      real, pointer, dimension(:,:,:) :: FSWUNAN, FSCUNAN

      real, pointer, dimension(:,:,:) ::    FSWD,    FSCD
      real, pointer, dimension(:,:,:) ::  FSWDNA,  FSCDNA

      real, pointer, dimension(:,:,:) ::     FSWBAND
      real, pointer, dimension(:,:,:) ::    FSWBANDN
      real, pointer, dimension(:,:,:) ::   FSWBANDNA
      real, pointer, dimension(:,:,:) ::  FSWBANDNAN

      real, pointer, dimension(:,:,:) ::  DRBANDN, DFBANDN
      real, pointer, dimension(:,:,:) ::  DROBIO , DFOBIO

      real, pointer, dimension(:,:  ) ::  DFUVR,  DRUVR
      real, pointer, dimension(:,:  ) ::  DFPAR,  DRPAR
      real, pointer, dimension(:,:  ) ::  DFNIR,  DRNIR
      real, pointer, dimension(:,:  ) :: DRNUVR, DRNPAR, DRNNIR
      real, pointer, dimension(:,:  ) ::  SLRTP,    RSR,  RSC
      real, pointer, dimension(:,:  ) ::  SLRSF,   RSCS,  RSRS, SLRSFC
      real, pointer, dimension(:,:  ) ::  SLRSFNA, SLRSFCNA
      real, pointer, dimension(:,:  ) ::  SLRSUF,  SLRSUFNA
      real, pointer, dimension(:,:  ) ::  SLRSUFC, SLRSUFCNA
      real, pointer, dimension(:,:  ) ::  RSRNA,  RSCNA
      real, pointer, dimension(:,:  ) :: RSCSNA, RSRSNA
      real, pointer, dimension(:,:  ) ::    OSR, OSRCLR
      real, pointer, dimension(:,:  ) ::  OSRNA, OSRCNA
      real, pointer, dimension(:,:  ) :: DFUVRN, DRUVRN
      real, pointer, dimension(:,:  ) :: DFPARN, DRPARN
      real, pointer, dimension(:,:  ) :: DFNIRN, DRNIRN
      real, pointer, dimension(:,:  ) :: ALBEDO
      real, pointer, dimension(:,:  ) :: COSZ, MCOSZ

      real, pointer, dimension(:,:,:)   :: FCLD,CLIN,RH
      real, pointer, dimension(:,:,:)   :: DP, PL, PLL, AERO, T, Q
      real, pointer, dimension(:,:,:)   :: RRL,RRI,RRR,RRS
      real, pointer, dimension(:,:,:)   :: RQL,RQI,RQR,RQS
      real, pointer, dimension(:,:,:,:) :: TAUCLD, HYDROMETS, REFF
      real, pointer, dimension(:,:,:)   :: TAUI,TAUW,TAUR,TAUS

      ! for efficiency
      real, allocatable, dimension(:,:) :: aCLDL,aCLDM,aCLDH
      real, allocatable, dimension(:,:) :: aTAUL,aTAUM,aTAUH
      real, allocatable, dimension(:,:) :: aCLDT,aTAUT

      real, dimension(LM  ) :: DUM1D
      real, dimension(LM,4) :: DUM2D

      real, pointer, dimension(:,:)   :: TDUST,TSALT,TSO4,TBC,TOC
      real, pointer, dimension(:,:)   :: CLDH,CLDM,CLDL,CLDT, &
                                         TAUH,TAUM,TAUL,TAUX,TAUT, &
                                         COTH,COTM,COTL,COTT, &
                                         CLDTMP,CLDPRS
      real, pointer, dimension(:,:)   :: COTDH,COTDM,COTDL,COTDT, &
                                         COTNH,COTNM,COTNL,COTNT

#ifdef SOLAR_RADVAL
      ! super-layer RRTMG cloud fraction exports on heartbeat
      real, pointer, dimension(:,:)   :: CLDTTSWHB
      real, pointer, dimension(:,:)   :: CLDHISWHB
      real, pointer, dimension(:,:)   :: CLDMDSWHB
      real, pointer, dimension(:,:)   :: CLDLOSWHB

      ! locals supporting CLD??SWHB
      integer :: rpart, pncol, ncld
      real    :: plmid(LM), tlev(LM-1), cfac(LM)
      integer, allocatable, dimension(:)     :: icld, jcld
      real,    allocatable, dimension(:)     :: alat
      real,    allocatable, dimension(:,:)   :: zmid, play, cldfrac, ciwp, clwp
      logical, allocatable, dimension(:,:,:) :: cldymcl
      real,    allocatable, dimension(:,:,:) :: ciwpmcl, clwpmcl
      integer, allocatable, dimension(:,:)   :: clearCounts
#endif

      type (ESMF_FieldBundle)         :: BUNDLE
      type (ESMF_Field)               :: FIELD
      type (ESMF_TimeInterval)        :: DELT
      character(len=ESMF_MAXSTR)      :: NAME

      integer :: L, I, J, N, NA, IB
      real    :: TAUCRIT
      real    :: FAC
      integer :: idx
      integer, external :: GetAeroIndex

      type (ty_RRTMGP_state), pointer      :: rrtmgp_state => null()
      type (ty_RRTMGP_wrap)                :: wrap
      character (len=ESMF_MAXPATHLEN)      :: k_dist_file
      character (len=ESMF_MAXSTR)          :: error_msg
      type (ty_gas_optics_rrtmgp), pointer :: k_dist
      type (ty_gas_concs)                  :: gas_concs

      ! OBIO bands (start,finish) in [nm]
      real, parameter :: OBIO_bands_nm (2,NB_OBIO) = reshape([ &
          200.0,  300.0, &  ! 01
          300.0,  350.0, &  ! 02
          350.0,  362.5, &  ! 03
          362.5,  387.5, &  ! 04
          387.5,  412.5, &  ! 05
          412.5,  437.5, &  ! 06
          437.5,  462.5, &  ! 07
          462.5,  487.5, &  ! 08
          487.5,  512.5, &  ! 09
          512.5,  537.5, &  ! 10
          537.5,  562.5, &  ! 11
          562.5,  587.5, &  ! 12
          587.5,  612.5, &  ! 13
          612.5,  637.5, &  ! 14
          637.5,  662.5, &  ! 15
          662.5,  687.5, &  ! 16
          687.5,  700.0, &  ! 17
          700.0,  750.0, &  ! 18
          750.0,  800.0, &  ! 19
          800.0,  900.0, &  ! 20
          900.0, 1000.0, &  ! 21
         1000.0, 1100.0, &  ! 22
         1100.0, 1200.0, &  ! 23
         1200.0, 1300.0, &  ! 24
         1300.0, 1400.0, &  ! 25
         1400.0, 1500.0, &  ! 26
         1500.0, 1600.0, &  ! 27
         1600.0, 1700.0, &  ! 28
         1700.0, 1800.0, &  ! 29
         1800.0, 2000.0, &  ! 30
         2000.0, 2400.0, &  ! 31
         2400.0, 3400.0, &  ! 32
         3400.0, 4000.0  &  ! 33
      ], shape(OBIO_bands_nm))

      ! We will do the conversion from SOLAR to OBIO bands in wavenumber [cm^-1],
      ! since photon energy \propto wavenumber, where wavenumber \def 1 / wavelength.
      ! The wavenumber bounds must clearly be swapped from the wavelength bounds
      ! because wavenumber is the inverse of wavelength.
      real, parameter :: OBIO_bands_wavenum (2,NB_OBIO) = 1.e7 / OBIO_bands_nm(2:1:-1,:)

      ! CHOU bands (start,finish) in [nm]
      real, parameter :: CHOU_bands_nm (2,NB_CHOU) = reshape([ &
          225.0,  285.0, &  ! 01
          285.0,  300.0, &  ! sub-band 2b (see note below)
          300.0,  325.0, &  ! 03
          325.0,  400.0, &  ! 04
          400.0,  690.0, &  ! 05
          690.0, 1220.0, &  ! 06
         1220.0, 2270.0, &  ! 07
         2270.0, 3850.0  &  ! 08
      ], shape(CHOU_bands_nm))

      ! Note: Chou-Suarez "Band 2" also includes sub-band 2a=(175,225) [nm].
      ! But almost no radiation in this extreme UV sub-band reaches the surface,
      ! so we treat all the band "2" as actually all sub-band 2b=(285,300) [nm],
      ! as in the CHOU_bands_nm array above.

      ! loaded later
      real :: SOLAR_bands_wavenum (2,NUM_BANDS_SOLAR)  ! [cm^-1]
      integer :: SOLAR_band_number_in_wvn_order (NUM_BANDS_SOLAR)

      real :: swvn1, swvn2, owvn1, owvn2, sfrac
      integer :: iseg, ibbeg, ibend, jb, kb, kb_start, kb_used_last
      logical :: sfirst, ofirst

      Iam  = trim(COMP_NAME)//"SolarUpdateExport"

      call ESMF_ClockGet(CLOCK, TIMESTEP=DELT, currTIME=CURRENTTIME, __RC__)

      call MAPL_SunGetInsolation(LONS, LATS, &
              ORBIT, ZTH, SLR, &
              INTV  = DELT,    &
              CLOCK = CLOCK,   &
              TIME = SUNFLAG,  &
              ZTHN = ZTHN,     &
              __RC__ )

      ZTH = max(ZTH,0.0)
      SLR = SLR * SC

      where(ZTH>0.0)
         SLN = (SLR/ZTH)
      else where
         SLN = 0.0
      end where

      call MAPL_GetPointer(INTERNAL, FSWN,       'FSWN',       __RC__)
      call MAPL_GetPointer(INTERNAL, FSCN,       'FSCN',       __RC__)
      call MAPL_GetPointer(INTERNAL, FSWNAN,     'FSWNAN',     __RC__)
      call MAPL_GetPointer(INTERNAL, FSCNAN,     'FSCNAN',     __RC__)
      call MAPL_GetPointer(INTERNAL, FSWUN,      'FSWUN',      __RC__)
      call MAPL_GetPointer(INTERNAL, FSCUN,      'FSCUN',      __RC__)
      call MAPL_GetPointer(INTERNAL, FSWUNAN,    'FSWUNAN',    __RC__)
      call MAPL_GetPointer(INTERNAL, FSCUNAN,    'FSCUNAN',    __RC__)
      call MAPL_GetPointer(INTERNAL, FSWBANDN,   'FSWBANDN',   __RC__)
      call MAPL_GetPointer(INTERNAL, FSWBANDNAN, 'FSWBANDNAN', __RC__)

      call MAPL_GetPointer(INTERNAL, DRUVRN,     'DRUVRN',     __RC__)
      call MAPL_GetPointer(INTERNAL, DFUVRN,     'DFUVRN',     __RC__)
      call MAPL_GetPointer(INTERNAL, DRPARN,     'DRPARN',     __RC__)
      call MAPL_GetPointer(INTERNAL, DFPARN,     'DFPARN',     __RC__)
      call MAPL_GetPointer(INTERNAL, DRNIRN,     'DRNIRN',     __RC__)
      call MAPL_GetPointer(INTERNAL, DFNIRN,     'DFNIRN',     __RC__)

      call MAPL_GetPointer(EXPORT  , FSW,        'FSW',        __RC__)
      call MAPL_GetPointer(EXPORT  , FSC,        'FSC',        __RC__)
      call MAPL_GetPointer(EXPORT  , FSWNA,      'FSWNA',      __RC__)
      call MAPL_GetPointer(EXPORT  , FSCNA,      'FSCNA',      __RC__)
      call MAPL_GetPointer(EXPORT  , FSWD,       'FSWD',       __RC__)
      call MAPL_GetPointer(EXPORT  , FSCD,       'FSCD',       __RC__)
      call MAPL_GetPointer(EXPORT  , FSWDNA,     'FSWDNA',     __RC__)
      call MAPL_GetPointer(EXPORT  , FSCDNA,     'FSCDNA',     __RC__)
      call MAPL_GetPointer(EXPORT  , FSWU,       'FSWU',       __RC__)
      call MAPL_GetPointer(EXPORT  , FSCU,       'FSCU',       __RC__)
      call MAPL_GetPointer(EXPORT  , FSWUNA,     'FSWUNA',     __RC__)
      call MAPL_GetPointer(EXPORT  , FSCUNA,     'FSCUNA',     __RC__)
      call MAPL_GetPointer(EXPORT  , FSWBAND,    'FSWBAND',    __RC__)
      call MAPL_GetPointer(EXPORT  , FSWBANDNA,  'FSWBANDNA',  __RC__)

      call MAPL_GetPointer(EXPORT  , DRUVR,      'DRUVR',      __RC__)
      call MAPL_GetPointer(EXPORT  , DFUVR,      'DFUVR',      __RC__)
      call MAPL_GetPointer(EXPORT  , DRPAR,      'DRPAR',      __RC__)
      call MAPL_GetPointer(EXPORT  , DFPAR,      'DFPAR',      __RC__)
      call MAPL_GetPointer(EXPORT  , DRNIR,      'DRNIR',      __RC__)
      call MAPL_GetPointer(EXPORT  , DFNIR,      'DFNIR',      __RC__)
      call MAPL_GetPointer(EXPORT  , RSR,        'RSR',        __RC__)
      call MAPL_GetPointer(EXPORT  , RSC,        'RSC',        __RC__)
      call MAPL_GetPointer(EXPORT  , RSRNA,      'RSRNA',      __RC__)
      call MAPL_GetPointer(EXPORT  , RSCNA,      'RSCNA',      __RC__)
      call MAPL_GetPointer(EXPORT  , SLRTP,      'SLRTP',      __RC__)
      call MAPL_GetPointer(EXPORT  , RSCS,       'RSCS',       __RC__)
      call MAPL_GetPointer(EXPORT  , RSRS,       'RSRS',       __RC__)
      call MAPL_GetPointer(EXPORT  , RSCSNA,     'RSCSNA',     __RC__)
      call MAPL_GetPointer(EXPORT  , RSRSNA,     'RSRSNA',     __RC__)
      call MAPL_GetPointer(EXPORT  , SLRSF,      'SLRSF',      __RC__)
      call MAPL_GetPointer(EXPORT  , SLRSFC,     'SLRSFC',     __RC__)
      call MAPL_GetPointer(EXPORT  , SLRSFNA,    'SLRSFNA',    __RC__)
      call MAPL_GetPointer(EXPORT  , SLRSFCNA,   'SLRSFCNA',   __RC__)
      call MAPL_GetPointer(EXPORT  , SLRSUF,     'SLRSUF',     __RC__)
      call MAPL_GetPointer(EXPORT  , SLRSUFC,    'SLRSUFC',    __RC__)
      call MAPL_GetPointer(EXPORT  , SLRSUFNA,   'SLRSUFNA',   __RC__)
      call MAPL_GetPointer(EXPORT  , SLRSUFCNA,  'SLRSUFCNA',  __RC__)
      call MAPL_GetPointer(EXPORT  , OSR,        'OSR',        __RC__)
      call MAPL_GetPointer(EXPORT  , OSRCLR,     'OSRCLR',     __RC__)
      call MAPL_GetPointer(EXPORT  , OSRNA,      'OSRNA',      __RC__)
      call MAPL_GetPointer(EXPORT  , OSRCNA,     'OSRCNA',     __RC__)
      call MAPL_GetPointer(EXPORT  , ALBEDO,     'ALBEDO',     __RC__)
      call MAPL_GetPointer(EXPORT  , COSZ,       'COSZ',       __RC__)
      call MAPL_GetPointer(EXPORT  , MCOSZ,      'MCOSZ',      __RC__)
      call MAPL_GetPointer(EXPORT  , DRNUVR,     'DRNUVR',     __RC__)
      call MAPL_GetPointer(EXPORT  , DRNPAR,     'DRNPAR',     __RC__)
      call MAPL_GetPointer(EXPORT  , DRNNIR,     'DRNNIR',     __RC__)

      call MAPL_GetPointer(IMPORT  , CLIN,       'FCLD',       __RC__)
      call MAPL_GetPointer(IMPORT  , PLL,        'PLE',        __RC__)
      call MAPL_GetPointer(IMPORT  , RRI,        'RI',         __RC__)
      call MAPL_GetPointer(IMPORT  , RRL,        'RL',         __RC__)
      call MAPL_GetPointer(IMPORT  , RRR,        'RR',         __RC__)
      call MAPL_GetPointer(IMPORT  , RRS,        'RS',         __RC__)
      call MAPL_GetPointer(IMPORT  , RQI,        'QI',         __RC__)
      call MAPL_GetPointer(IMPORT  , RQL,        'QL',         __RC__)
      call MAPL_GetPointer(IMPORT  , RQR,        'QR',         __RC__)
      call MAPL_GetPointer(IMPORT  , RQS,        'QS',         __RC__)
      call MAPL_GetPointer(IMPORT  , T,          'T',          __RC__)
      call MAPL_GetPointer(IMPORT  , Q,          'QV',         __RC__)

      call MAPL_GetPointer(EXPORT  , FCLD,       'FCLD',       __RC__)
      call MAPL_GetPointer(EXPORT  , TAUI,       'TAUCLI',     __RC__)
      call MAPL_GetPointer(EXPORT  , TAUW,       'TAUCLW',     __RC__)
      call MAPL_GetPointer(EXPORT  , TAUR,       'TAUCLR',     __RC__)
      call MAPL_GetPointer(EXPORT  , TAUS,       'TAUCLS',     __RC__)
      call MAPL_GetPointer(EXPORT  , CLDL,       'CLDLO',      __RC__)
      call MAPL_GetPointer(EXPORT  , CLDM,       'CLDMD',      __RC__)
      call MAPL_GetPointer(EXPORT  , CLDH,       'CLDHI',      __RC__)
      call MAPL_GetPointer(EXPORT  , CLDT,       'CLDTT',      __RC__)
      call MAPL_GetPointer(EXPORT  , TAUL,       'TAULO',      __RC__)
      call MAPL_GetPointer(EXPORT  , TAUM,       'TAUMD',      __RC__)
      call MAPL_GetPointer(EXPORT  , TAUH,       'TAUHI',      __RC__)
      call MAPL_GetPointer(EXPORT  , TAUT,       'TAUTT',      __RC__)
      call MAPL_GetPointer(EXPORT  , TAUX,       'TAUTX',      __RC__)
      call MAPL_GetPointer(EXPORT  , COTL,       'COTLO',      __RC__)
      call MAPL_GetPointer(EXPORT  , COTM,       'COTMD',      __RC__)
      call MAPL_GetPointer(EXPORT  , COTH,       'COTHI',      __RC__)
      call MAPL_GetPointer(EXPORT  , COTT,       'COTTT',      __RC__)
      call MAPL_GetPointer(EXPORT  , CLDTMP,     'CLDTMP',     __RC__)
      call MAPL_GetPointer(EXPORT  , CLDPRS,     'CLDPRS',     __RC__)
      call MAPL_GetPointer(EXPORT  , COTDL,      'COTDENLO',   __RC__)
      call MAPL_GetPointer(EXPORT  , COTDM,      'COTDENMD',   __RC__)
      call MAPL_GetPointer(EXPORT  , COTDH,      'COTDENHI',   __RC__)
      call MAPL_GetPointer(EXPORT  , COTDT,      'COTDENTT',   __RC__)
      call MAPL_GetPointer(EXPORT  , COTNL,      'COTNUMLO',   __RC__)
      call MAPL_GetPointer(EXPORT  , COTNM,      'COTNUMMD',   __RC__)
      call MAPL_GetPointer(EXPORT  , COTNH,      'COTNUMHI',   __RC__)
      call MAPL_GetPointer(EXPORT  , COTNT,      'COTNUMTT',   __RC__)

#ifdef SOLAR_RADVAL
      call MAPL_GetPointer(EXPORT  , CLDLOSWHB,  'CLDLOSWHB',  __RC__)
      call MAPL_GetPointer(EXPORT  , CLDMDSWHB,  'CLDMDSWHB',  __RC__)
      call MAPL_GetPointer(EXPORT  , CLDHISWHB,  'CLDHISWHB',  __RC__)
      call MAPL_GetPointer(EXPORT  , CLDTTSWHB,  'CLDTTSWHB',  __RC__)
#endif

      if (SOLAR_TO_OBIO) then
         call MAPL_GetPointer(INTERNAL, DRBANDN, 'DRBANDN',    __RC__)
         call MAPL_GetPointer(INTERNAL, DFBANDN, 'DFBANDN',    __RC__)
         call MAPL_GetPointer(EXPORT  , DROBIO,  'DROBIO',     __RC__)
         call MAPL_GetPointer(EXPORT  , DFOBIO,  'DFOBIO',     __RC__)
      end if

      if (associated(FCLD)) FCLD = CLIN

      if (associated(CLDH) .or. associated(CLDT) .or. &
          associated(TAUX) .or. associated(COTT) .or. &
          associated(COTDH) .or. associated(COTNH) .or. &
          associated(COTDT) .or. associated(COTNT)) &
      then
         allocate(aCLDH(IM,JM),__STAT__)
         aCLDH = 0.
         do l=1,LCLDMH-1
            aCLDH = max(aCLDH,CLIN(:,:,L))
         end do
         if (associated(CLDH)) CLDH = aCLDH
         if (associated(COTDH)) COTDH = aCLDH
      end if

      if (associated(CLDM) .or. associated(CLDT) .or. &
          associated(TAUX) .or. associated(COTT) .or. &
          associated(COTDM) .or. associated(COTNM) .or. &
          associated(COTDT) .or. associated(COTNT)) &
      then
         allocate(aCLDM(IM,JM),__STAT__)
         aCLDM = 0.
         do l=LCLDMH,LCLDLM-1
            aCLDM = max(aCLDM,CLIN(:,:,L))
         end do
         if (associated(CLDM)) CLDM = aCLDM
         if (associated(COTDM)) COTDM = aCLDM
      end if

      if (associated(CLDL) .or. associated(CLDT) .or. &
          associated(TAUX) .or. associated(COTT) .or. &
          associated(COTDL) .or. associated(COTNL) .or. &
          associated(COTDT) .or. associated(COTNT)) &
      then
         allocate(aCLDL(IM,JM),__STAT__)
         aCLDL = 0.
         do l=LCLDLM,LM
            aCLDL = max(aCLDL,CLIN(:,:,L))
         end do
         if (associated(CLDL)) CLDL = aCLDL
         if (associated(COTDL)) COTDL = aCLDL
      end if

      if (associated(CLDT) .or. &
          associated(TAUX) .or. associated(COTT) .or. &
          associated(COTDT) .or. associated(COTNT)) &
      then
         allocate(aCLDT(IM,JM),__STAT__)
         aCLDT = 1. - (1-aCLDH)*(1-aCLDM)*(1-aCLDL)
         if (associated(CLDT)) CLDT = aCLDT
         if (associated(COTDT)) COTDT = aCLDT
      end if

#ifdef SOLAR_RADVAL
      ! CLD??SWHB:
      ! Special heartbeat versions of RRTMG generated cloud fractions ...
      ! These are expensive because they require a call to the cloud generator,
      ! which normally is only done inside the IRRAD and SOLAR REFRESHes (and
      ! for the SOLAR case only on the sunlit portion of the globe). We provide
      ! these here as a means of validation, but they should not be regularly
      ! exported since they will slow down SOLAR considerably, and since the
      ! equivalent REFRESH-generated versions (without the "HB" suffix), which
      ! are updated only at the REFRESH frequency (~hourly), should be fine in
      ! most cases (especially for longer term averages).
      ! NB: Filled on all globe unlike the RESFRESH version CLD??SW.

      if (associated(CLDLOSWHB) .or. associated(CLDMDSWHB) .or. &
          associated(CLDHISWHB) .or. associated(CLDTTSWHB)) then

         ! default to clear columns which do not need subcolumn generation
         if (associated(CLDLOSWHB)) CLDLOSWHB = 0.
         if (associated(CLDMDSWHB)) CLDMDSWHB = 0.
         if (associated(CLDHISWHB)) CLDHISWHB = 0.
         if (associated(CLDTTSWHB)) CLDTTSWHB = 0.

         ! partition size pncol for cloudy columns to conserve memory & improve efficiency
         call MAPL_GetResource(MAPL,rpart,'RRTMGSW_PARTITION_SIZE:',DEFAULT=0,__RC__)
         if (rpart > 0) then
            pncol = rpart
         else
            pncol = 2
         end if

         ! space for partition:
         ! The partition stores up cloudy gridcolumns to process in batch
         allocate(icld   (          pncol),__STAT__)
         allocate(jcld   (          pncol),__STAT__)
         allocate(zmid   (LM,       pncol),__STAT__)
         allocate(alat   (          pncol),__STAT__)
         allocate(play   (LM,       pncol),__STAT__)
         allocate(cldfrac(LM,       pncol),__STAT__)
         allocate(ciwp   (LM,       pncol),__STAT__)
         allocate(clwp   (LM,       pncol),__STAT__)
         allocate(cldymcl(LM,ngptsw,pncol),__STAT__)
         allocate(ciwpmcl(LM,ngptsw,pncol),__STAT__)
         allocate(clwpmcl(LM,ngptsw,pncol),__STAT__)
         allocate(clearCounts(4,    pncol),__STAT__)

         ! start with empty partition
         ncld = 0

! after this test ... make sure DOY used consistently by refresh and update in model
! I guess its only now being used in update, but was its use in refresh really consistent?
! this may be a non-zero-diff bug fix later ... co2 by DOY = hb but DOY for RRTMG should be for REFRESH style time

         ! loop over domain
         do j = 1,JM
            do i = 1,IM

               ! load up cloudy columns to partition
               if (any(CLIN(i,j,:) > 0.)) then

                  ! cloudy column
                  ncld = ncld + 1
                  icld (ncld) = i
                  jcld (ncld) = j
                  alat (ncld) = LATS(i,j)

                  ! Note: unlike RRTMG we do not reverse the levels. This is a
                  ! technicality and will not alter the POPULATION stats of the
                  ! generation and saves time (see notes under cloud_subcol_gen).
                  ! If an exact replication of RRTMGSW is required, can reverse
                  ! vertical ordering here ... but an exact replication will
                  ! also require saving the exact cldfrac used by the REFRESH
                  ! into the internal state for use here as well.

                  plmid = 0.5 * (PLL(i,j,0:LM-1) + PLL(i,j,1:LM))
                  play   (:,ncld) = plmid / 100.  ! hPa
                  cldfrac(:,ncld) = CLIN(i,j,:)

                  ! cloud water paths converted from g/g to g/m^2
                  cfac = 1.02 * 100 * (PLL(i,j,1:LM)-PLL(i,j,0:LM-1))
                  ciwp(:,ncld) = cfac * RQI(i,j,:)
                  clwp(:,ncld) = cfac * RQL(i,j,:)

                  ! interior interface temperatures
                  ! * "0-based" but extema at 0 and LM not needed for zmid;
                  ! * RRTMG call code uses layer delP (DPR) but any multiple
                  ! of it, specifically cfac, is equivalent.

                  tlev = (T(i,j,1:LM-1) * cfac(2:LM) + T(i,j,2:LM) * cfac(1:LM-1)) &
                       / (                cfac(2:LM) +               cfac(1:LM-1))

                  ! Calculate the LAYER (mid-point) heights.
                  ! The interlayer distances are needed for the calculations
                  ! of inter-layer correlation for cloud overlapping. Only
                  ! *relative* distances matter, so wolog set zmid(LM) = 0.
                  zmid(LM,ncld) = 0.
                  do k = LM-1, 1, -1
                     ! dz ~ RT/g x dp/p by hydrostatic eqn and ideal gas eqn.
                     ! The jump from LAYER k+1 to k is centered on LEVEL k
                     !   since the LEVEL indices are zero-based
                     zmid(k,ncld) = zmid(k+1,ncld) + MAPL_RGAS * tlev(k) / MAPL_GRAV &
                                        * (plmid(k+1) - plmid(k)) / PLL(i,j,k)
                  end do

               end if ! cloudy column

               ! nothing to process yet?
               if (ncld == 0) cycle

               ! process the partition if its full or if its partially
               !   full but there are no more gridcolumns left.

               if (ncld == pncol .or. i == IM .and. j == JM) then

                  ! McICA subcolumn generation
                  call generate_stochastic_clouds( &
                     pncol, ncld, ngptsw, LM, &
                     zmid, alat, doy, &
                     play, cldfrac, ciwp, clwp, 1.e-20, &
                     cldymcl, ciwpmcl, clwpmcl, &
                     seed_order=[4,3,2,1])

                  ! for super-layer cloud fractions
                  call clearCounts_threeBand( &
                     pncol, ncld, ngptsw, LM, LCLDLM, LCLDMH, cldymcl, &
                     clearCounts)

                  ! convert super-layer clearCounts to cloud fractions
                  if (associated(CLDTTSWHB)) then
                     do n = 1,ncld
                        CLDTTSWHB(icld(n),jcld(n)) = 1. - clearCounts(1,n)/float(ngptsw)
                     end do
                  end if
                  if (associated(CLDHISWHB)) then
                     do n = 1,ncld
                        CLDHISWHB(icld(n),jcld(n)) = 1. - clearCounts(2,n)/float(ngptsw)
                     end do
                  end if
                  if (associated(CLDMDSWHB)) then
                     do n = 1,ncld
                        CLDMDSWHB(icld(n),jcld(n)) = 1. - clearCounts(3,n)/float(ngptsw)
                     end do
                  end if
                  if (associated(CLDLOSWHB)) then
                     do n = 1,ncld
                        CLDLOSWHB(icld(n),jcld(n)) = 1. - clearCounts(4,n)/float(ngptsw)
                     end do
                  end if

                  ! restart partition
                  ncld = 0

               end if  ! process partition

            end do  ! i
         end do  ! j

         ! clean up
         deallocate(icld,jcld,__STAT__)
         deallocate(zmid,alat,play,cldfrac,ciwp,clwp,__STAT__)
         deallocate(cldymcl,ciwpmcl,clwpmcl,__STAT__)
         deallocate(clearCounts,__STAT__)

      end if  ! CLD??SWHB
#endif

      if (associated(TAUI) .or. associated(TAUW) .or. associated(TAUR) .or. associated(TAUS).or. &
          associated(TAUL) .or. associated(TAUM) .or. associated(TAUH) .or. &
          associated(COTL) .or. associated(COTM) .or. associated(COTH) .or. &
          associated(TAUT) .or. associated(TAUX) .or. associated(COTT) .or. &
          associated(COTNL) .or. associated(COTNM) .or. associated(COTNH) .or. associated(COTNT) .or. &
          associated(CLDTMP) .or. associated(CLDPRS)) &
      then

         allocate(   TAUCLD(IM,JM,LM,4), __STAT__)
         allocate(HYDROMETS(IM,JM,LM,4), __STAT__)
         allocate(     REFF(IM,JM,LM,4), __STAT__)
         allocate(       DP(IM,JM,LM  ), __STAT__)

         DP = PLL(:,:,1:LM)-PLL(:,:,0:LM-1)

         ! In REFF, HYDROMETS, and TAUCLD, the final index is as follows:
         !       1  Cloud Ice
         !       2  Cloud Liquid
         !       3  Falling Liquid (Rain)
         !       4  Falling Ice (Rain)

         REFF(:,:,:,1) = RRI * 1.e6  ! REFF must be in microns
         REFF(:,:,:,2) = RRL * 1.e6
         REFF(:,:,:,3) = RRR * 1.e6
         REFF(:,:,:,4) = RRS * 1.e6

         HYDROMETS(:,:,:,1) = RQI
         HYDROMETS(:,:,:,2) = RQL
         HYDROMETS(:,:,:,3) = RQR
         HYDROMETS(:,:,:,4) = RQS

         TAUCLD = 0.

         ! Due to the generic use of this routine, it currently works on one column at a time,
         ! thus the need for the array sections below.

         ! NOTE: Dummy arrays are passed into outputs 1 and 3 because these are currently only
         !       used in sorad.F90.

         DO I = 1,IM
            DO J = 1,JM
               CALL GETVISTAU( &
                  LM,ZTH(I,J),DP(I,J,:),&
                  CLIN(I,J,:),REFF(I,J,:,:),HYDROMETS(I,J,:,:),&
                  LCLDMH,LCLDLM,&
                  DUM2D(:,:),TAUCLD(I,J,:,:),DUM1D(:))
            END DO
         END DO

         if (associated(TAUI)) TAUI = TAUCLD(:,:,:,1)
         if (associated(TAUW)) TAUW = TAUCLD(:,:,:,2)
         if (associated(TAUR)) TAUR = TAUCLD(:,:,:,3)
         if (associated(TAUS)) TAUS = TAUCLD(:,:,:,4)

         ! use the total hydrometor optical thickness for the general opticl thicknesses below
         TAUCLD(:,:,:,1) = TAUCLD(:,:,:,1) + TAUCLD(:,:,:,2) + TAUCLD(:,:,:,3) + TAUCLD(:,:,:,4)

         ! TAU[HML] are correct because GETVISTAU produces in-cloud optical thicknesses for
         ! 'effective clouds' extended-out and diluted to the maximum cloud fraction in each
         ! pressure super-layers [LMH].

         if (associated(TAUH) .or. associated(COTH) .or. associated(COTNH) .or. &
             associated(TAUT) .or. associated(TAUX) .or. associated(COTT) .or. associated(COTNT)) &
         then
            allocate(aTAUH(IM,JM),__STAT__)
            aTAUH = 0.
            do l=1,LCLDMH-1
               aTAUH = aTAUH + TAUCLD(:,:,L,1)
            end do
            if (associated(TAUH)) TAUH = aTAUH
            if (associated(COTH)) then
              COTH = MAPL_UNDEF
              where (aCLDH > 0.) COTH = aTAUH
            end if
            if (associated(COTNH)) COTNH = aCLDH * aTAUH
         end if

         if (associated(TAUM) .or. associated(COTM) .or. associated(COTNM) .or. &
             associated(TAUT) .or. associated(TAUX) .or. associated(COTT) .or. associated(COTNT)) &
         then
            allocate(aTAUM(IM,JM),__STAT__)
            aTAUM = 0.
            do l=LCLDMH,LCLDLM-1
               aTAUM = aTAUM + TAUCLD(:,:,L,1)
            end do
            if (associated(TAUM)) TAUM = aTAUM
            if (associated(COTM)) then
              COTM = MAPL_UNDEF
              where (aCLDM > 0.) COTM = aTAUM
            end if
            if (associated(COTNM)) COTNM = aCLDM * aTAUM
         end if

         if (associated(TAUL) .or. associated(COTL) .or. associated(COTNL) .or. &
             associated(TAUT) .or. associated(TAUX) .or. associated(COTT) .or. associated(COTNT)) &
         then
            allocate(aTAUL(IM,JM),__STAT__)
            aTAUL = 0.
            do l=LCLDLM,LM
               aTAUL = aTAUL + TAUCLD(:,:,L,1)
            end do
            if (associated(TAUL)) TAUL = aTAUL
            if (associated(COTL)) then
              COTL = MAPL_UNDEF
              where (aCLDL > 0.) COTL = aTAUL
            end if
            if (associated(COTNL)) COTNL = aCLDL * aTAUL
         end if

         ! TAUT however is broken because the three super-layers are randomly overlapped
         ! and with different effective cloud fractions. It has been broken but used for
         ! a long time. It should be considered deprecated. TAUX below is an improved
         ! version.

         if (associated(TAUT)) TAUT = aTAUH + aTAUM + aTAUL

         ! As noted above, one cannot simply add TAUL, TAUM and TAUH to get a column
         ! in-cloud optical thickness, because the actual column value depends on the
         ! overlap of these bands. This overlap is here assumed random. We can express
         ! the approximate column in-cloud optical thickness in terms of the sum over
         ! the 2**3 - 1 combinations with some cloud in at least one of the 3 bands,
         ! each with their respective fractions. For random overlap, CLDL*CLDM*CLDH of
         ! the gridcolumn would have a column TAU of TAUL+TAUM+TAUH, CLDL*CLDM*(1-CLDH)
         ! would have a column TAU of TAUL+TAUM, etc. Then, for an in-cloud column TAU,
         ! the sum of the 7 must be normalized by the random column cloud fraction
         !    CLDT = 1 – (1-CLDL)*(1-CLDM)*(1-CLDH).
         ! Not surprisingly this gives
         !    TAUX = (TAUL*CLDL + TAUM*CLDM + TAUH*CLDH) / CLDT,
         ! because we assume we can linearly average optical thickness among the comb-
         ! inations. This assumption is questionable, since cloud radiative properties
         ! are non-linear in optical thickness. This is why TAUX is approximate. But
         ! its the best we SIMPLY can do.

         if (associated(TAUX) .or. associated(COTT) .or. associated(COTNT)) then
            allocate(aTAUT(IM,JM),__STAT__)
            aTAUT = 0.
            where (aCLDT > 0.) aTAUT = (aTAUL*aCLDL + aTAUM*aCLDM + aTAUH*aCLDH) / aCLDT
            if (associated(TAUX)) TAUX = aTAUT
            if (associated(COTT)) then
              COTT = MAPL_UNDEF 
              where (aCLDT > 0.) COTT = aTAUT
            end if
            if (associated(COTNT)) COTNT = aCLDT * aTAUT
         end if

         if (allocated(aTAUH)) deallocate(aTAUH,__STAT__)
         if (allocated(aTAUM)) deallocate(aTAUM,__STAT__)
         if (allocated(aTAUL)) deallocate(aTAUL,__STAT__)
         if (allocated(aTAUT)) deallocate(aTAUT,__STAT__)

         if (associated(CLDTMP) .or. associated(CLDPRS)) then
            call MAPL_GetResource(MAPL,TAUCRIT,'TAUCRIT:',DEFAULT=0.10,__RC__)

            if (associated(CLDTMP)) CLDTMP = MAPL_UNDEF
            if (associated(CLDPRS)) CLDPRS = MAPL_UNDEF

            do L=LM,1,-1
               if (associated(CLDTMP)) then
                  where (TAUCLD(:,:,L,1) > TAUCRIT) CLDTMP = T(:,:,L)
               end if
               if (associated(CLDPRS)) then
                  where (TAUCLD(:,:,L,1) > TAUCRIT) CLDPRS = PLL(:,:,L-1)
               end if
            end do
         end if

         deallocate(TAUCLD   )
         deallocate(HYDROMETS)
         deallocate(REFF     )
         deallocate(DP       )

      end if

      if (allocated(aCLDH)) deallocate(aCLDH,__STAT__)
      if (allocated(aCLDM)) deallocate(aCLDM,__STAT__)
      if (allocated(aCLDL)) deallocate(aCLDL,__STAT__)
      if (allocated(aCLDT)) deallocate(aCLDT,__STAT__)

! Fill Albedos
!-------------

      FAC = 1.

! Visible/UV diffuse

      call MAPL_GetPointer(EXPORT, ALBEXP, 'ALBVF', __RC__)
      if(associated(ALBEXP)) then
         call MAPL_GetPointer(IMPORT, ALBIMP, 'ALBVF', __RC__)
         where(SLR>0)
            ALBEXP = ALBIMP * FAC
         elsewhere
            ALBEXP = MAPL_UNDEF
         end where
      end if

! Visible/UV direct

      call MAPL_GetPointer(EXPORT, ALBEXP, 'ALBVR', __RC__)
      if(associated(ALBEXP)) then
         call MAPL_GetPointer(IMPORT, ALBIMP, 'ALBVR', __RC__)
         where(SLR>0)
            ALBEXP = ALBIMP * FAC
         elsewhere
            ALBEXP = MAPL_UNDEF
         end where
      end if

! NIR diffuse

      call MAPL_GetPointer(EXPORT, ALBEXP, 'ALBNF', __RC__)
      if(associated(ALBEXP)) then
         call MAPL_GetPointer(IMPORT, ALBIMP, 'ALBNF', __RC__)
         where(SLR>0)
            ALBEXP = ALBIMP * FAC
         elsewhere
            ALBEXP = MAPL_UNDEF
         end where
      end if

! NIR direct

      call MAPL_GetPointer(EXPORT, ALBEXP, 'ALBNR', __RC__)
      if(associated(ALBEXP)) then
         call MAPL_GetPointer(IMPORT, ALBIMP, 'ALBNR', __RC__)
         where(SLR>0)
            ALBEXP = ALBIMP * FAC
         elsewhere
            ALBEXP = MAPL_UNDEF
         end where
      end if

! Total surface albedo
!---------------------

      ALB = DRUVRN+DFUVRN+DRPARN+DFPARN+DRNIRN+DFNIRN
      where(SLR>0.0 .and. ALB>0.0)
         ALB = min( max(1.0 - FSWN(:,:,LM)/ALB,.01), 0.9 )
      elsewhere
         ALB = MAPL_UNDEF
      end where

      if(associated(ALBEDO)) ALBEDO = ALB

! Fill incident fluxes
!---------------------

      if(associated( SLRTP))  SLRTP =        SLR
      if(associated( DRUVR))  DRUVR = DRUVRN*SLR
      if(associated( DFUVR))  DFUVR = DFUVRN*SLR
      if(associated( DRPAR))  DRPAR = DRPARN*SLR
      if(associated( DFPAR))  DFPAR = DFPARN*SLR
      if(associated( DRNIR))  DRNIR = DRNIRN*SLR
      if(associated( DFNIR))  DFNIR = DFNIRN*SLR
      if(associated(DRNUVR)) DRNUVR = DRUVRN*SLN
      if(associated(DRNPAR)) DRNPAR = DRPARN*SLN
      if(associated(DRNNIR)) DRNNIR = DRNIRN*SLN

      if(associated( SLRSF))  SLRSF = (DRUVRN+DFUVRN+DRPARN+DFPARN+DRNIRN+DFNIRN)*SLR

      if(associated(SLRSFC)) then
         where(ALB /= MAPL_UNDEF)
            SLRSFC = (FSCN(:,:,LM)*SLR)/(1.-ALB)
         elsewhere
            SLRSFC = 0.0
         end where
      end if

      if(associated(SLRSFNA)) then
         where(ALB /= MAPL_UNDEF)
            SLRSFNA = (FSWNAN(:,:,LM)*SLR)/(1.-ALB)
         elsewhere
            SLRSFNA = 0.0
         end where
      end if

      if(associated(SLRSFCNA)) then
         where(ALB /= MAPL_UNDEF)
            SLRSFCNA = (FSCNAN(:,:,LM)*SLR)/(1.-ALB)
         elsewhere
            SLRSFCNA = 0.0
         end where
      end if

      if(associated(   SLRSUF)) SLRSUF  = ALB*(DRUVRN+DFUVRN+DRPARN+DFPARN+DRNIRN+DFNIRN)*SLR

      if(associated(  SLRSUFC)) then
        where(ALB /= MAPL_UNDEF)
           SLRSUFC = ALB*(FSCN(:,:,LM)/(1.-ALB))*SLR
        elsewhere
           SLRSUFC = 0.0
        end where
      end if

      if(associated( SLRSUFNA)) then
        where(ALB /= MAPL_UNDEF)
           SLRSUFNA = ALB*(FSWNAN(:,:,LM)/(1.-ALB))*SLR
        elsewhere
           SLRSUFNA = 0.0
        end where
      end if

      if(associated(SLRSUFCNA)) then
        where(ALB /= MAPL_UNDEF)
           SLRSUFCNA = ALB*(FSCNAN(:,:,LM)/(1.-ALB))*SLR
        elsewhere
           SLRSUFCNA = 0.0
        end where
      end if

! Fill 3D FLuxes
!---------------

      DO L=0,LM
         ! Fill Export Net Fluxes from Internal
         ! ------------------------------------
         if(associated(FSW  )) FSW  (:,:,L) = FSWN  (:,:,L)*SLR
         if(associated(FSC  )) FSC  (:,:,L) = FSCN  (:,:,L)*SLR
         if(associated(FSWNA)) FSWNA(:,:,L) = FSWNAN(:,:,L)*SLR
         if(associated(FSCNA)) FSCNA(:,:,L) = FSCNAN(:,:,L)*SLR

         ! Fill Export Up Fluxes from Internal
         ! -----------------------------------
         if(associated(FSWU  )) FSWU  (:,:,L) = FSWUN  (:,:,L)*SLR
         if(associated(FSCU  )) FSCU  (:,:,L) = FSCUN  (:,:,L)*SLR
         if(associated(FSWUNA)) FSWUNA(:,:,L) = FSWUNAN(:,:,L)*SLR
         if(associated(FSCUNA)) FSCUNA(:,:,L) = FSCUNAN(:,:,L)*SLR

         ! Fill Export Down Fluxes from (Net + Up) Internal
         ! ------------------------------------------------
         if(associated(FSWD  )) FSWD  (:,:,L) = (FSWN  (:,:,L) + FSWUN  (:,:,L))*SLR
         if(associated(FSCD  )) FSCD  (:,:,L) = (FSCN  (:,:,L) + FSCUN  (:,:,L))*SLR
         if(associated(FSWDNA)) FSWDNA(:,:,L) = (FSWNAN(:,:,L) + FSWUNAN(:,:,L))*SLR
         if(associated(FSCDNA)) FSCDNA(:,:,L) = (FSCNAN(:,:,L) + FSCUNAN(:,:,L))*SLR
      end do

! Fill 3D per-band Fluxes
! -----------------------

      do IB = 1, NUM_BANDS_SOLAR
         if(associated(FSWBAND  )) FSWBAND  (:,:,IB) = FSWBANDN  (:,:,IB)*SLR
         if(associated(FSWBANDNA)) FSWBANDNA(:,:,IB) = FSWBANDNAN(:,:,IB)*SLR
      end do

! Fill 2D Fluxes
!---------------

      if(associated(   RSR))    RSR =       FSWN(:,:, 0) *SLR
      if(associated(   RSC))    RSC =       FSCN(:,:, 0) *SLR
      if(associated( RSRNA))  RSRNA =     FSWNAN(:,:, 0) *SLR
      if(associated( RSCNA))  RSCNA =     FSCNAN(:,:, 0) *SLR
      if(associated(  RSRS))   RSRS =       FSWN(:,:,LM) *SLR
      if(associated(  RSCS))   RSCS =       FSCN(:,:,LM) *SLR
      if(associated(RSRSNA)) RSRSNA =     FSWNAN(:,:,LM) *SLR
      if(associated(RSCSNA)) RSCSNA =     FSCNAN(:,:,LM) *SLR
      if(associated(   OSR))    OSR = (1.-  FSWN(:,:, 0))*SLR
      if(associated(OSRCLR)) OSRCLR = (1.-  FSCN(:,:, 0))*SLR
      if(associated( OSRNA))  OSRNA = (1.-FSWNAN(:,:, 0))*SLR
      if(associated(OSRCNA)) OSRCNA = (1.-FSCNAN(:,:, 0))*SLR

! SOLAR TO OBIO conversion ...
! Done in wavenum [cm^-1] space for reasons detailed under OBIO_bands_wavenum declaration
! ---------------------------------------------------------------------------------------

      if (SOLAR_TO_OBIO) then
         if (associated(DROBIO) .or. associated(DFOBIO)) then

            ! first load SOLAR_bands_wavenum and specify ordering
            !    that makes it monotonically increasing ...

            if (USE_RRTMGP) then

! helper for testing RRTMGP error status on return;
! allows line number reporting cf. original call method
#define TEST_(A) error_msg = A; if (trim(error_msg)/="") then; _FAIL("RRTMGP Error: "//trim(error_msg)); endif

               ! access RRTMGP internal state from the GC
               call ESMF_UserCompGetInternalState(GC, 'RRTMGP_state', wrap, status)
               VERIFY_(status)
               rrtmgp_state => wrap%ptr

               ! initialize k-distribution if not already done
               ! remember: its possible to have UPDATE_FIRST
               if (.not. rrtmgp_state%initialized) then
                  call MAPL_GetResource( &
                     MAPL, k_dist_file, "RRTMGP_DATA_SW:", &
                     DEFAULT='rrtmgp-data-sw.nc', __RC__)
                  ! gas_concs needed only to access required gas names
                  error_msg = gas_concs%init([character(3) :: &
                     'h2o','co2','o3','n2o','co','ch4','o2','n2'])
                  TEST_(error_msg)
                  call load_and_init( &
                     rrtmgp_state%k_dist, trim(k_dist_file), gas_concs)
                  if (.not. rrtmgp_state%k_dist%source_is_external()) then
                     TEST_('RRTMGP-SW: does not seem to be SW')
                  endif
                  rrtmgp_state%initialized = .true.
               endif

               ! access by shorter name
               k_dist => rrtmgp_state%k_dist

               ! load RRTMGP bands (2,NB_RRTMGP) [cm^-1]
               SOLAR_bands_wavenum = k_dist%get_band_lims_wavenumber()

               ! RRTMGP bands are already ordered in increasing wavenumber
               SOLAR_band_number_in_wvn_order = [(i, i=1,NUM_BANDS_SOLAR)]

#undef TEST_

            elseif (USE_RRTMG) then

               ! note: RRTMG bands are [wavenum1(jpb1:jpb2),wavenum2(jpb1:jpb2)] in [cm^-1]
               ! The index jpb1:jpb2 (16:29) is over the 14 bands. Band 14 is OUT of order.

               _ASSERT(jpb2-jpb1+1     == 14, 'RRTMG band inconsistency!')
               _ASSERT(NUM_BANDS_SOLAR == 14, 'wrong number of RRTMG bands!')

               ! load RRTMG bands (2,NB_RRTMG) [cm^-1]
               SOLAR_bands_wavenum (1,:) = wavenum1(jpb1:jpb2)
               SOLAR_bands_wavenum (2,:) = wavenum2(jpb1:jpb2)

               ! RRTMG band 14 comes before band 1 in increasing wavenumber
               SOLAR_band_number_in_wvn_order(1) = 14
               SOLAR_band_number_in_wvn_order(2:14) = [(i, i=1,13)]

            else if (USE_CHOU) then

               ! Source CHOU_bands_nm (2,NB_CHOU) is ordered in increasing waveLENGTH.
               ! See CHOU_bands_nm declaration for note about band 2a/b.

               ! Load Chou-Suarez bands (2,NB_CHOU) [cm^-1].
               ! Flip waveLENGTH bounds to waveNUMBER bounds.
               SOLAR_bands_wavenum (1,:) = 1.e7 / CHOU_bands_nm(2,:)
               SOLAR_bands_wavenum (2,:) = 1.e7 / CHOU_bands_nm(1,:)

               ! Specify sorting of Chou-Suarez bands in increasing waveNUMBER
               SOLAR_band_number_in_wvn_order = [(i, i=NUM_BANDS_SOLAR,1,-1)]

            endif

            ! zero accumulators
            if (associated(DROBIO)) DROBIO = 0.
            if (associated(DFOBIO)) DFOBIO = 0.

            ! We attempt to do the conversion efficiently by taking into account the band
            ! ordering in increasing wavenumber and knowing there are no gaps or overlaps
            ! in either SOLAR or OBIO bands.
            sfirst = .true.
            ofirst = .true.
            kb_start = NB_OBIO
            do jb = 1,NUM_BANDS_SOLAR
               ib = SOLAR_band_number_in_wvn_order(jb)

               ! SOLAR band (swvn1,swvn2)
               ! Note: check SOLAR band continuity before updating swvn2
               swvn1 = SOLAR_bands_wavenum(1,ib)
               if (.not.sfirst) then
                  _ASSERT(swvn1 == swvn2, 'SOLAR bands not complete and unique!')
               end if
               swvn2 = SOLAR_bands_wavenum(2,ib)
               sfirst = .false.

               ! Now find which OBIO bands the SOLAR band contributes to. OBIO bands
               ! are in increasing wavelength, decreasing wavenumber, so this loop is
               ! in increasing wavenumber.
               do kb = kb_start,1,-1

                  ! OBIO band (owvn1,owvn2)
                  ! Note: check OBIO band continuity before updating owvn2
                  owvn1 = OBIO_bands_wavenum(1,kb)
                  if (.not.ofirst) then
                     if (kb .ne. kb_used_last) then
                        _ASSERT(owvn1 == owvn2, 'OBIO bands not complete and unique!')
                     end if
                  end if
                  owvn2 = OBIO_bands_wavenum(2,kb)
                  kb_used_last = kb
                  ofirst = .false.

                  ! if we exit the OBIO band loop for any reason, begin
                  ! processing the next SOLAR band into the same OBIO band
                  kb_start = kb

                  ! OBIO band has wavenumbers all higher than current SOLAR
                  ! band, so immediately move on to the next SOLAR band
                  if (owvn1 >= swvn2) exit

                  ! OBIO band has wavenumbers all lower than current SOLAR
                  ! band, so skip it and move on to the next OBIO band
                  if (owvn2 <= swvn1) cycle

                  ! now there is some overlap between the SOLAR and OBIO bands ...

                  ! accumulate assuming constant energy spread across each SOLAR band, sfrac in (0,1]
                  sfrac = (min(swvn2,owvn2)-max(swvn1,owvn1)) / (swvn2 - swvn1)
                  if (associated(DROBIO)) DROBIO (:,:,kb) = DROBIO (:,:,kb) + DRBANDN (:,:,ib) * sfrac
                  if (associated(DFOBIO)) DFOBIO (:,:,kb) = DFOBIO (:,:,kb) + DFBANDN (:,:,ib) * sfrac

                  ! OBIO band has some wavenumbers higher than current
                  ! SOLAR band, so move on to the next SOLAR band
                  if (owvn2 > swvn2) exit

               end do  ! kb (OBIO band)
            end do  ! jb (SOLAR band)

            ! unnormalise to W/m2
            do kb = 1,NB_OBIO
               if (associated(DROBIO)) DROBIO (:,:,kb) = DROBIO (:,:,kb) * SLR
               if (associated(DFOBIO)) DFOBIO (:,:,kb) = DFOBIO (:,:,kb) * SLR
            end do

         end if
      end if

! Solar zenith angles: mean for time step and at end of time step.
!  Note SLR should not be used after this point!!
!-----------------------------------------------------------------

      if(associated( MCOSZ))  MCOSZ = ZTH
      if(associated(  COSZ))   COSZ = ZTHN

      RETURN_(ESMF_SUCCESS)
    end subroutine UPDATE_EXPORT


  end subroutine RUN

  ! Pack masked locations into buffer
  subroutine PackIt (Packed, UnPacked, MSK, Pdim, Udim, LM)
    integer, intent(IN   ) :: Pdim, Udim(2), LM
    real,    intent(INOUT) ::   Packed(Pdim,*)
    real,    intent(IN   ) :: UnPacked(Udim(1),Udim(2),*)
    logical, intent(IN   ) :: MSK(Udim(1),Udim(2))

    integer :: I, J, L, M

    do L = 1,LM
      M = 1
      do J = 1,Udim(2)
        do I = 1,Udim(1)
          if (MSK(I,J)) then
            Packed(M,L) = UnPacked(I,J,L)
            M = M+1
          end if
        end do
      end do
    end do

  end subroutine PackIt

  ! Unpack masked locations from buffer
  subroutine UnPackIt(Packed, UnPacked, MSK, Pdim, Udim, LM, DEFAULT)
    integer, intent(IN   ) :: Pdim, Udim(2), LM
    real,    intent(IN   ) ::   Packed(Pdim,*)
    real,    intent(INOUT) :: UnPacked(Udim(1),Udim(2),*)
    logical, intent(IN   ) :: MSK(Udim(1),Udim(2))
    real, optional, intent(IN) :: DEFAULT

    integer :: I, J, L, M

    do L = 1,LM
      M = 1
      do J = 1,Udim(2)
        do I = 1,Udim(1)
          if (MSK(I,J)) then
            Unpacked(I,J,L) = Packed(M,L)
            M = M+1
          elseif (PRESENT(DEFAULT)) then
            UnPacked(I,J,L) = DEFAULT
          end if
        end do
      end do
    end do

  end subroutine UnPackIt

  ! Decide which radiation to use for thermodynamics state evolution.
  ! RRTMGP dominates RRTMG dominates Chou-Suarez.
  ! Chou-Suarez is the default if nothing else asked for in Resource file.
  !----------------------------------------------------------------------

  subroutine choose_solar_scheme (MAPL, &
    USE_RRTMGP, USE_RRTMG, USE_CHOU, &
    RC)

    type (MAPL_MetaComp), pointer, intent(in) :: MAPL
    logical, intent(out) :: USE_RRTMGP, USE_RRTMG, USE_CHOU
    integer, optional, intent(out) :: RC  ! return code

    real :: RFLAG
    integer :: STATUS

    USE_RRTMGP = .false.
    USE_RRTMG  = .false.
    USE_CHOU   = .false.
    call MAPL_GetResource (MAPL, RFLAG, LABEL='USE_RRTMGP_SORAD:', DEFAULT=0., __RC__)
    USE_RRTMGP = RFLAG /= 0.
    if (.not. USE_RRTMGP) then
      call MAPL_GetResource (MAPL, RFLAG, LABEL='USE_RRTMG_SORAD:', DEFAULT=0., __RC__)
      USE_RRTMG = RFLAG /= 0.
      USE_CHOU  = .not.USE_RRTMG
    end if

    _RETURN(_SUCCESS)
  end subroutine choose_solar_scheme


  subroutine choose_irrad_scheme (MAPL, &
    USE_RRTMGP, USE_RRTMG, USE_CHOU, &
    RC)

    type (MAPL_MetaComp), pointer, intent(in) :: MAPL
    logical, intent(out) :: USE_RRTMGP, USE_RRTMG, USE_CHOU
    integer, optional, intent(out) :: RC  ! return code

    real :: RFLAG
    integer :: STATUS

    USE_RRTMGP = .false.
    USE_RRTMG  = .false.
    USE_CHOU   = .false.
    call MAPL_GetResource (MAPL, RFLAG, LABEL='USE_RRTMGP_IRRAD:', DEFAULT=0., __RC__)
    USE_RRTMGP = RFLAG /= 0.
    if (.not. USE_RRTMGP) then
      call MAPL_GetResource (MAPL, RFLAG, LABEL='USE_RRTMG_IRRAD:', DEFAULT=0., __RC__)
      USE_RRTMG = RFLAG /= 0.
      USE_CHOU  = .not.USE_RRTMG
    end if

    _RETURN(_SUCCESS)
  end subroutine choose_irrad_scheme

end module GEOS_SolarGridCompMod
