# - Try to find Ipopt
# Once done this will define
# Ipopt_FOUND           - system has SuiteSparse
# Ipopt_INCLUDE_DIRS    - the SuiteSparse include directories
# Ipopt_LIBRARIES       - Link these to use SuiteSparse
#

if (NOT Ipopt_FOUND)
  find_path(Ipopt_INCLUDE_DIR IpTNLP.hpp
            PATHS
            "${IPOPT_HOME}/include/coin"
            "${Ipopt_HOME}/include/coin"
            "$ENV{IPOPT_HOME}/include/coin"
            "/usr/include/coin"
            NO_DEFAULT_PATH)
  find_library(Ipopt_LIBRARY ipopt
               PATHS
               "${IPOPT_HOME}/lib"
               "${Ipopt_HOME}/lib"
               "$ENV{IPOPT_HOME}/lib"
               "/usr/lib"
               NO_DEFAULT_PATH)


  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(Ipopt DEFAULT_MSG Ipopt_LIBRARY Ipopt_INCLUDE_DIR)

  if (Ipopt_FOUND)
    set(Ipopt_INCLUDE_DIRS "${Ipopt_INCLUDE_DIR}")
    set(Ipopt_LIBRARIES "${Ipopt_LIBRARY}")

    add_library(Ipopt::Ipopt INTERFACE IMPORTED)
    set_target_properties(Ipopt::Ipopt PROPERTIES INTERFACE_LINK_LIBRARIES ${Ipopt_LIBRARIES})
    set_target_properties(Ipopt::Ipopt PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${Ipopt_INCLUDE_DIRS})
  endif ()

  mark_as_advanced(Ipopt_INCLUDE_DIRS Ipopt_LIBRARIES)
endif ()

