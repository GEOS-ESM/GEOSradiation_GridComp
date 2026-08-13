add_subdirectory(FakeGOCART_GridComp)

file(STRINGS "cases.txt" TEST_CASES)

include(esma_regression_tests)
esma_add_regression_tests(
  NAME irrad
  DATA_PATH GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSradiation_GridComp
  TEST_CASES ${TEST_CASES}
)
