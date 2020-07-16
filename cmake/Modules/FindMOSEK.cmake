#
# Try to find MOSEK
# Once done this will define
#
# MOSEK_FOUND           - system has MOSEK
# MOSEK_INCLUDE_DIRS    - the MOSEK include directories
# MOSEK_LIBRARIES       - Link these to use MOSEK
#

#if already found skip the search
IF (NOT MOSEK_FOUND)
  SET(SEARCH_PATHS
      /usr/local/
      /usr/
      "${MOSEK_LIBRARY_DIR}"
      "${MOSEKROOT}"
      "${MOSEK_HOME}"
      "$ENV{MOSEKROOT}"
      "$ENV{MOSEK_HOME}"
      )

  FIND_PATH(MOSEK_INCLUDE_DIR mosek.h
            PATHS ${SEARCH_PATHS}
            PATH_SUFFIXES h)

  FIND_LIBRARY(MOSEK_LIBRARY NAMES mosek64
               PATHS ${SEARCH_PATHS}
               PATH_SUFFIXES bin)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(MOSEK
                                    FOUND_VAR MOSEK_FOUND
                                    REQUIRED_VARS
                                    MOSEK_LIBRARY
                                    MOSEK_INCLUDE_DIR
                                    )

  mark_as_advanced(MOSEK_INCLUDE_DIR MOSEK_LIBRARIES)

  IF (MOSEK_FOUND)
    set(MOSEK_LIBRARIES ${MOSEK_LIBRARY})
    set(MOSEK_INCLUDE_DIRS ${MOSEK_INCLUDE_DIR})

    add_library(MOSEK::MOSEK INTERFACE IMPORTED)
    set_target_properties(MOSEK::MOSEK PROPERTIES INTERFACE_LINK_LIBRARIES ${MOSEK_LIBRARIES})
    set_target_properties(MOSEK::MOSEK PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${MOSEK_INCLUDE_DIRS})
  ENDIF (MOSEK_FOUND)



  mark_as_advanced(MOSEK_INCLUDE_DIR MOSEK_LIBRARIES)
endif ()
