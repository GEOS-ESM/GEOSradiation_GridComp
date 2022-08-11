esma_set_this()

set (alldirs
  GEOS_SolarGridComp
  GEOS_IrradGridComp
  GEOS_SatsimGridComp
  )

esma_add_subdirectories(
  Shared
  )

esma_add_library (${this}
  SRCS GEOS_RadiationGridComp.F90
  SUBCOMPONENTS ${alldirs}
  DEPENDENCIES MAPL GEOS_Shared esmf GEOS_RadiationShared)

