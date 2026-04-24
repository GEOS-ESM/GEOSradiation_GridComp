esma_set_this()

set (srcs
  soradconstants.F90
  sorad.F90
  GEOS_SolarGridComp.F90
  )

esma_add_library (${this}
  SRCS ${srcs}
  SUBDIRS RRTMG
  DEPENDENCIES MAPL GEOS_Shared GEOS_RadiationShared RRTMG_SW RRTMG_SW_mods RRTMGP)

if(MKL_FOUND)
   target_compile_definitions(${this} PRIVATE HAVE_MKL)
endif ()

option(ENABLE_SOLAR_RADVAL "extra solar diagnostics mainly for developer validation" OFF)
if(ENABLE_SOLAR_RADVAL)
   target_compile_definitions(${this} PRIVATE SOLAR_RADVAL)
endif ()
