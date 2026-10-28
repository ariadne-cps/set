find_library(ARIADNE_SET_LIBRARY NAMES ariadne-set)

find_package(Threads REQUIRED)

find_path(ARIADNE_SET_INCLUDE_DIR ariadne-set.hpp PATH_SUFFIXES ariadne-set)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(AriadneSet DEFAULT_MSG
  ARIADNE_SET_LIBRARY
  ARIADNE_SET_INCLUDE_DIR
)

if(AriadneSet_FOUND)
  get_filename_component(ARIADNE_SET_INCLUDE_PARENT_DIR ${ARIADNE_SET_INCLUDE_DIR} DIRECTORY)
  set(ARIADNE_SET_INCLUDE_DIRS
    ${ARIADNE_SET_INCLUDE_PARENT_DIR}
    ${ARIADNE_SET_INCLUDE_DIR}
  )
  if(NOT TARGET AriadneSet::ariadne-set)
    add_library(AriadneSet::ariadne-set UNKNOWN IMPORTED)
    set_target_properties(AriadneSet::ariadne-set PROPERTIES
      IMPORTED_LOCATION "${ARIADNE_SET_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${ARIADNE_SET_INCLUDE_DIRS}"
      INTERFACE_LINK_LIBRARIES "Threads::Threads"
    )
  endif()
endif()
