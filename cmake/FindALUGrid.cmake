# - Try to find ALUGrid

find_path(ALUGRID_INCLUDE_DIR alugrid_parallel.h HINTS ${ALUGRID_ROOT}/include)

find_library(ALUGRID_LIBRARY alugrid HINTS ${ALUGRID_ROOT}/lib)

set(ALUGRID_LIBRARIES ${ALUGRID_LIBRARY})
set(ALUGRID_INCLUDE_DIRS ${ALUGRID_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)


find_package_handle_standard_args(ALUGRID DEFAULT_MSG ALUGRID_LIBRARIES ALUGRID_INCLUDE_DIRS)

