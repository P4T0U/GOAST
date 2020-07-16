cmake_minimum_required(VERSION 3.5)

# project name
project(GOAST VERSION 0.1)
set(CMAKE_CXX_STANDARD 11)

################################################################################
# general settings
################################################################################

# where to find the modules
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)

################################################################################
# Library CMake include
################################################################################
include(cmake/GOAST.cmake)

###########################################
# Use MKL (via Eigen) in the examples etc
###########################################
option(GOAST_WITH_MKL "Use Intel's Math Kernel Library" OFF)
if (GOAST_WITH_MKL)
  # target_link_options(GOAST INTERFACE -m64 -Wl,--no-as-needed)
  # interface link flags not supported in CMake 3.5, hence
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-as-needed")

  target_compile_options(GOAST INTERFACE -m64)
  target_include_directories(GOAST INTERFACE $ENV{MKLROOT}/include)

  target_link_directories(GOAST INTERFACE $ENV{MKLROOT}/lib/intel64)
  target_link_libraries(GOAST INTERFACE mkl_intel_lp64 mkl_gnu_thread mkl_core gomp pthread m dl)

  target_compile_definitions(GOAST INTERFACE GOAST_WITH_MKL)
endif (GOAST_WITH_MKL)

################################################################################
# source directories
################################################################################
option(GOAST_BUILD_TESTS "Build GOAST tests" ON)
option(GOAST_BUILD_EXAMPLES "Build examples of GOAST" ON)
option(GOAST_BUILD_TOOLS "Build tools included with GOAST" ON)

if (GOAST_BUILD_TESTS)
  add_subdirectory(tests)
endif ()

if (GOAST_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif ()

if (GOAST_BUILD_TOOLS)
  add_subdirectory(tools)
endif ()


################################################################################
# install directories
################################################################################
include(GNUInstallDirs)


install(TARGETS GOAST GOAST_all
        EXPORT GOASTTargets
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        )

install(EXPORT GOASTTargets
        FILE GOASTTargets.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/GOAST
        )

include(CMakePackageConfigHelpers)
write_basic_package_version_file("GOASTConfigVersion.cmake"
                                 VERSION ${GOAST_VERSION}
                                 COMPATIBILITY SameMajorVersion)

configure_file(cmake/GOASTConfig.cmake.in GOASTConfig.cmake @ONLY)

install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/GOASTConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/GOASTConfigVersion.cmake"
        "cmake/Modules/FindSuiteSparse.cmake"
        "cmake/Modules/FindADOLC.cmake"
        "cmake/Modules/FindEigen3.cmake"
        "cmake/Modules/FindIpopt.cmake"
        "cmake/Modules/FindMOSEK.cmake"
        "cmake/Modules/FindOpenMesh.cmake"
        "cmake/Modules/Findtrlib.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/GOAST)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/goast DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})