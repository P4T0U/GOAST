project(examples)

add_executable(multires_example multires_example.cpp)
add_executable(interpolationShells geodesic_interpolation_shells.cpp)
add_executable(extrapolationShells geodesic_extrapolation_shells.cpp)
add_executable(dissipationModes dissipationModes.cpp)
add_executable(decimation decimation.cpp)
add_executable(elasticMean elasticMean.cpp)
add_executable(parallelTransport parallelTransport.cpp)
add_executable(deformationTransfer deformationTransfer.cpp)
add_executable(nonlinearSpline nonlinearSpline.cpp)
add_executable(sectionalCurvature sectionalCurvature.cpp)
add_executable(registration registration.cpp)
add_executable(DirichletEnergy DirichletEnergy.cpp)
add_executable(parametricMCM parametricMCM.cpp)

# link libraries
target_link_libraries(multires_example PRIVATE GOAST::GOAST)

target_link_libraries(interpolationShells PRIVATE GOAST::GOAST)

target_link_libraries(extrapolationShells PRIVATE GOAST::GOAST)

target_link_libraries(dissipationModes PRIVATE GOAST::GOAST)

target_link_libraries(decimation PRIVATE GOAST::GOAST)

target_link_libraries(elasticMean PRIVATE GOAST::GOAST)

target_link_libraries(parallelTransport PRIVATE GOAST::GOAST)

target_link_libraries(deformationTransfer PRIVATE GOAST::GOAST)

target_link_libraries(nonlinearSpline PRIVATE GOAST::GOAST)

target_link_libraries(sectionalCurvature PRIVATE GOAST::GOAST)

target_link_libraries(registration PRIVATE GOAST::GOAST)

target_link_libraries(DirichletEnergy PRIVATE GOAST::GOAST)

target_link_libraries(parametricMCM PRIVATE GOAST::GOAST)

# Phasefield example currently requires Ipopt
if (GOAST_WITH_IPOPT)
    add_executable(Phasefield Phasefield.cpp)
    target_link_libraries(Phasefield PRIVATE GOAST)
    target_link_libraries(Phasefield PRIVATE GOAST::All)
endif ()

# NRIC Examples
add_executable(NRICGeodesics NRICGeodesics.cpp)
target_link_libraries(NRICGeodesics PRIVATE GOAST::GOAST)

add_executable(NRICConstrictions NRICConstrictions.cpp)
target_link_libraries(NRICConstrictions PRIVATE GOAST::GOAST)

add_executable(NRICFolding NRICFolding.cpp)
target_link_libraries(NRICFolding PRIVATE GOAST::GOAST)

add_executable(NRICMean NRICMean.cpp)
target_link_libraries(NRICMean PRIVATE GOAST::GOAST)

add_executable(NRICIsometricMode NRICIsometricMode.cpp)
target_link_libraries(NRICIsometricMode PRIVATE GOAST::All)

add_executable(nonlinearProjection nonlinearProjection.cpp)
target_link_libraries(nonlinearProjection PRIVATE GOAST::GOAST)