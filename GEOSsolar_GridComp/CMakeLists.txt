esma_set_this()

set (srcs
  soradconstants.F90
  sorad.F90
  GEOS_SolarGridComp.F90
  )

esma_add_library (${this}
  SRCS ${srcs}
  SUBDIRS RRTMG
  DEPENDENCIES MAPL GEOS_Shared GEOS_RadiationShared RRTMG_SW RRTMG_SW_mods RRTMGP fvdycore)

if(MKL_FOUND)
   target_compile_definitions(${this} PRIVATE HAVE_MKL)
endif ()

option(ENABLE_SOLAR_RADVAL "extra solar diagnostics mainly for developer validation" OFF)
if(ENABLE_SOLAR_RADVAL)
   target_compile_definitions(${this} PRIVATE SOLAR_RADVAL)
endif ()

# Install GEOS-MLT ML radiation runtime files.
set(MLRAD_RUNTIME_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}/@GEOSmlrad_GridComp
)

install(
    FILES
        ${MLRAD_RUNTIME_DIR}/python/geos_mlrad_driver.py
    DESTINATION
        lib/Python
)

install(
    FILES
        ${MLRAD_RUNTIME_DIR}/models/mlrad_colfilm1d_qrs.pth
        ${MLRAD_RUNTIME_DIR}/models/mlrad_colfilm1d_qrl.pth
        ${MLRAD_RUNTIME_DIR}/models/mlrad_colfilm1d_qjoule.pth
    DESTINATION
        share/GEOSmlrad/models
)

install(
    FILES
        ${MLRAD_RUNTIME_DIR}/norms/NormalizationStats_GEOSgrid_2000-2010.npz
    DESTINATION
        share/GEOSmlrad/norms
)
