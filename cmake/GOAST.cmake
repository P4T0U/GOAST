cmake_minimum_required(VERSION 3.5)

################################################################################
# Essential paths
################################################################################
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/Modules)

################################################################################
# Mandatory libraries: Eigen, OpenMesh, and SuiteSparse
################################################################################
if (NOT TARGET Eigen3::Eigen)
  find_package(Eigen3 REQUIRED)
endif ()
if (NOT TARGET OpenMeshCore)
  find_package(OpenMesh REQUIRED)
endif ()
if (NOT TARGET SuiteSparse::suitesparseconfig)
  find_package(SuiteSparse REQUIRED)
endif ()

################################################################################
# GOAST targets
################################################################################
# Basic GOAST
add_library(GOAST INTERFACE)
add_library(GOAST::GOAST ALIAS GOAST)
set_property(TARGET GOAST PROPERTY EXPORT_NAME GOAST::GOAST)


target_include_directories(GOAST INTERFACE
                           "$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../include>"
                           "$<INSTALL_INTERFACE:include>")

target_link_libraries(GOAST INTERFACE
                      Eigen3::Eigen
                      SuiteSparse::suitesparseconfig
                      OpenMeshCore
                      OpenMeshTools)

# Target for GOAST with all optional dependencies included
add_library(GOAST_all INTERFACE)
add_library(GOAST::All ALIAS GOAST_all)
set_property(TARGET GOAST_all PROPERTY EXPORT_NAME GOAST::All)
target_link_libraries(GOAST_all INTERFACE GOAST) # GOAST with all extension includes regular GOAST

################################################################################
# OpenMP
################################################################################
option(GOAST_WITH_OPENMP "Enable use of OpenMP within GOAST" ON)

if (GOAST_WITH_OPENMP)
  find_package(OpenMP)

  #for cmake version less than 3.9 there is no find module, so we need to define the target ourselves
  if (NOT TARGET OpenMP::OpenMP_CXX)
    if (${CMAKE_CXX_COMPILER_ID} STREQUAL "AppleClang" OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
      message(WARNING "OpenMP detection might now work with AppleClang, Intel")
    endif ()
    find_package(Threads REQUIRED)
    add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
    set_property(TARGET OpenMP::OpenMP_CXX
                 PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
    set_property(TARGET OpenMP::OpenMP_CXX
                 PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_CXX_FLAGS} Threads::Threads)
  endif ()

  target_link_libraries(GOAST INTERFACE OpenMP::OpenMP_CXX)
  target_compile_definitions(GOAST INTERFACE GOAST_WITH_OPENMP)
endif (GOAST_WITH_OPENMP)

################################################################################
# Optional libraries
################################################################################
option(GOAST_WITH_IPOPT "Enable Ipopt interface withing GOAST" OFF)
option(GOAST_WITH_ADOLC "Use the ADOL-C autodiff package" OFF)
option(GOAST_WITH_VTK "Use the VTK library within GOAST" OFF)
option(GOAST_WITH_TRLIB "Use the Trlib trust-region subproblem solver package within GOAST" OFF)
option(GOAST_WITH_MOSEK "Enable MOSEK interface within GOAST" OFF)

# todo: move to imported targets to make it packagable
###########################################
# IPOPT
###########################################
if (GOAST_WITH_IPOPT)
  find_package(Ipopt REQUIRED)

  target_link_libraries(GOAST_all INTERFACE Ipopt::Ipopt)
  target_compile_definitions(GOAST_all INTERFACE GOAST_WITH_IPOPT)
  target_compile_definitions(GOAST_all INTERFACE HAVE_CSTDDEF) # wrong config under Debian workaround
endif ()

###########################################
# ADOL-C
###########################################
if (GOAST_WITH_ADOLC)
  find_package(ADOLC REQUIRED)

  target_link_libraries(GOAST_all INTERFACE ADOLC::ADOLC)
  target_compile_definitions(GOAST_all INTERFACE GOAST_WITH_ADOLC)
endif ()

###########################################
# VTK
###########################################
if (GOAST_WITH_VTK)
  find_package(VTK REQUIRED)
  if (VTK_VERSION VERSION_LESS "8.90.0")
    include(${VTK_USE_FILE})
  endif ()

  target_include_directories(GOAST_all SYSTEM INTERFACE ${VTK_INCLUDE_DIRS})
  target_link_libraries(GOAST_all INTERFACE ${VTK_LIBRARIES})

  target_compile_definitions(GOAST_all INTERFACE GOAST_WITH_VTK)
endif ()

###########################################
# trlib
###########################################
if (GOAST_WITH_TRLIB)
  find_package(trlib REQUIRED)

  target_link_libraries(GOAST_all INTERFACE trlib::trlib)
  target_compile_definitions(GOAST_all INTERFACE GOAST_WITH_TRLIB)
endif ()

###########################################
# MOSEK
###########################################
if (GOAST_WITH_MOSEK)
  find_package(MOSEK MODULE REQUIRED)

  target_link_libraries(GOAST_all INTERFACE MOSEK::MOSEK)
  target_compile_definitions(GOAST_all INTERFACE GOAST_WITH_MOSEK)
endif ()
