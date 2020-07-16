project(tests)

add_executable(RosenbrockTest RosenbrockTest.cpp)
target_link_libraries(RosenbrockTest PRIVATE GOAST)

add_executable(NRICIntegrabilityTest NRICIntegrabilityTest.cpp)
target_link_libraries(NRICIntegrabilityTest PRIVATE GOAST::All)

add_executable(NRICReconstructionTest NRICReconstructionTest.cpp)
target_link_libraries(NRICReconstructionTest PRIVATE GOAST::All)

if (GOAST_WITH_ADOLC)
  add_executable(adolcTest adolcTest.cpp)
  target_link_libraries(adolcTest PRIVATE GOAST::All)
endif ()

add_executable(bosqpTest bosqpTest.cpp)
target_link_libraries(bosqpTest PRIVATE GOAST::All)

add_executable(firstTest firstTest.cpp)
target_link_libraries(firstTest PRIVATE GOAST::GOAST)

if (GOAST_WITH_IPOPT)
  add_executable(ipoptTest ipoptTest.cpp)
  target_link_libraries(ipoptTest PRIVATE GOAST::All)
endif ()

if (GOAST_WITH_MOSEK)
  add_executable(MOSEKTest MOSEKTest.cpp)
  target_link_libraries(MOSEKTest PRIVATE GOAST::All)
endif ()

add_executable(trncgTest trncgTest.cpp)
target_link_libraries(trncgTest PRIVATE GOAST::All)