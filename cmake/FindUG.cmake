# - Try to find UG
# Once done this will define
#  LIBXML2_FOUND - System has LibXml2
#  LIBXML2_INCLUDE_DIRS - The LibXml2 include directories
#  LIBXML2_LIBRARIES - The libraries needed to use LibXml2
#  LIBXML2_DEFINITIONS - Compiler switches required for using LibXml2

#find_package(PkgConfig)
#pkg_check_modules(UG ugS2)
#set(LIBXML2_DEFINITIONS ${PC_LIBXML_CFLAGS_OTHER})
set (UG_ROOT /home/mrueckl/opt/ug)

find_path(UG_INCLUDE_DIR ug/ugenv.h HINTS ${UG_ROOT}/include)

find_library(UG_S2_LIBRARY ugS2 HINTS ${UG_ROOT}/lib)
find_library(UG_S3_LIBRARY ugS3 HINTS ${UG_ROOT}/lib)
find_library(UG_devS_LIBRARY devS HINTS ${UG_ROOT}/lib)

set(UG_LIBRARIES ${UG_S2_LIBRARY} ${UG_S3_LIBRARY} ${UG_devS_LIBRARY})
set(UG_INCLUDE_DIRS ${UG_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBXML2_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(UG DEFAULT_MSG UG_S2_LIBRARY UG_S2_LIBRARY)
find_package_handle_standard_args(UG DEFAULT_MSG UG_S3_LIBRARY UG_S3_LIBRARY)
find_package_handle_standard_args(UG DEFAULT_MSG UG_devS_LIBRARY UG_devS_LIBRARY)

#mark_as_advanced(UG_INCLUDE_DIR UG_LIBRARY)
