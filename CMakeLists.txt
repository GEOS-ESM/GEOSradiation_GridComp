esma_set_this()

set (alldirs
  GEOS_RadiationShared
  # GEOSsolar_GridComp NOT ported to MAPL3 yet
  GEOSirrad_GridComp
  # GEOSsatsim_GridComp
)

esma_add_library (${this}
  SRCS GEOS_RadiationGridComp.F90
  SUBCOMPONENTS ${alldirs}
  DEPENDENCIES MAPL GEOS_Shared ESMF::ESMF
  TYPE SHARED
)

mapl_acg (
  ${this} Radiation_StateSpecs.rc
  IMPORT_SPECS EXPORT_SPECS GET_POINTERS DECLARE_POINTERS
)

add_subdirectory(regression)

