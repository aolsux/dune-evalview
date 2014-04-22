# - Try to find Metis
find_path(METIS_INCLUDE_DIR metis.h HINTS ${METIS_ROOT}/include )

find_library(METIS_LIBRARY metis HINTS ${METIS_ROOT}/lib )

set(METIS_LIBRARIES ${METIS_LIBRARY})
set(METIS_INCLUDE_DIRS ${METIS_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Metis DEFAULT_MSG METIS_INCLUDE_DIRS METIS_LIBRARIES)

