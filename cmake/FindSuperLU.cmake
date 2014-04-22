# - Try to find UG
# Once done this will define
#  LIBXML2_FOUND - System has LibXml2
#  LIBXML2_INCLUDE_DIRS - The LibXml2 include directories
#  LIBXML2_LIBRARIES - The libraries needed to use LibXml2
#  LIBXML2_DEFINITIONS - Compiler switches required for using LibXml2

#find_package(PkgConfig)
#pkg_check_modules(UG ugS2)
#set(LIBXML2_DEFINITIONS ${PC_LIBXML_CFLAGS_OTHER})

find_path(SuperLU_INCLUDE_DIR slu_cdefs.h HINTS /users/tsd/mrueckl/install/SuperLU_4.3/SRC)

find_library(SuperLU_LIBRARY superlu_4.3 HINTS /users/tsd/mrueckl/install/SuperLU_4.3/lib)

set(SuperLU_LIBRARIES ${SuperLU_LIBRARY} -lcblas -lf77blas -latlas)
set(SuperLU_INCLUDE_DIRS ${SuperLU_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBXML2_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(SuperLU DEFAULT_MSG SuperLU_LIBRARIES SuperLU_INCLUDE_DIR)
set(SuperLU_FOUND TRUE)
#mark_as_advanced(UG_INCLUDE_DIR UG_LIBRARY)
