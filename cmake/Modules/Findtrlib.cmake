#
# Try to find trlib
# Once done this will define
#  
# trlib_FOUND           - system has trlib
# trlib_INCLUDE_DIR     - the trlib include directory
# trlib_LIBRARIES         - Link this to use trlib
# trlib_LIBRARY_DIR     - directory where the libraries are included
#
#===========================================================================

cmake_minimum_required(VERSION 2.8.9)

#if already found skip the search
IF (NOT trlib_FOUND)
  SET(SEARCH_PATHS
      /usr/local/
      /usr/
      "${trlib_LIBRARY_DIR}"
      "${trlib_HOME}"
      "${TRLIB_HOME}"
      "$ENV{TRLIB_HOME}"
      )

  FIND_PATH(trlib_INCLUDE_DIR trlib.h
      PATHS ${SEARCH_PATHS}
      PATH_SUFFIXES include)

  FIND_LIBRARY(trlib_LIBRARIES NAMES trlib
      PATHS ${SEARCH_PATHS}
      PATH_SUFFIXES lib)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(trlib DEFAULT_MSG trlib_LIBRARIES trlib_INCLUDE_DIR)

  get_filename_component(_trlib_LIBRARY_DIR ${trlib_LIBRARIES} PATH)
  set(trlib_LIBRARY_DIR "${_trlib_LIBRARY_DIR}" CACHE PATH "The directory where the trlib library can be found. ")

  add_library(trlib::trlib INTERFACE IMPORTED)
  set_target_properties(trlib::trlib PROPERTIES INTERFACE_LINK_LIBRARIES ${trlib_LIBRARIES})
  set_target_properties(trlib::trlib PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${trlib_INCLUDE_DIR})

  mark_as_advanced(trlib_INCLUDE_DIR trlib_LIBRARIES trlib_LIBRARY_DIR)
endif ()
