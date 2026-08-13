esma_set_this()

esma_add_library (${this}
  SRCS FakeGOCART_GridComp.F90
  DEPENDENCIES MAPL ESMF::ESMF
  TYPE SHARED)
