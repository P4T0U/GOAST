# - Try to find ADOL-C
# Once done this will define
# ADOLC_FOUND           - system has SuiteSparse
# ADOLC_INCLUDE_DIRS    - the SuiteSparse include directories
# ADOLC_LIBRARIES       - Link these to use SuiteSparse
#

if (NOT ADOLC_FOUND)
  find_path(ADOLC_INCLUDE_DIR
            NAMES
            adolc/adtl.h
            PATHS
            "${ADOLC_HOME}/include"
            "$ENV{ADOLC_HOME}/include"
            "$ENV{ADOLCDIR}/include"
            ${INCLUDE_INSTALL_DIR}
            "/usr/include"
            NO_DEFAULT_PATH)

  find_library(ADOLC_LIBRARY
               NAMES
               adolc
               PATHS
               "${ADOLC_HOME}"
               "$ENV{ADOLC_HOME}"
               "$ENV{ADOLCDIR}"
               ${LIB_INSTALL_DIR}
               "/usr"
               PATH_SUFFIXES lib lib64
               NO_DEFAULT_PATH)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(ADOLC DEFAULT_MSG ADOLC_LIBRARY ADOLC_INCLUDE_DIR)

  if (ADOLC_FOUND)
    set(ADOLC_INCLUDE_DIRS "${ADOLC_INCLUDE_DIR}")
    set(ADOLC_LIBRARIES "${ADOLC_LIBRARY}")

    add_library(ADOLC::ADOLC INTERFACE IMPORTED)
    set_target_properties(ADOLC::ADOLC PROPERTIES INTERFACE_LINK_LIBRARIES ${ADOLC_LIBRARIES})
    set_target_properties(ADOLC::ADOLC PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${ADOLC_INCLUDE_DIRS})
  endif ()

  mark_as_advanced(ADOLC_INCLUDE_DIRS ADOLC_LIBRARIES)
endif ()