# - Try to find DUNE

find_path(DUNE_COMMON_INCLUDE_DIR         dune/common         HINTS ${DUNE_ROOT}/include ${DUNE_ROOT}/dune-common)
find_path(DUNE_GEOMETRY_INCLUDE_DIR       dune/geometry       HINTS ${DUNE_ROOT}/include ${DUNE_ROOT}/dune-geometry)
find_path(DUNE_GRID_INCLUDE_DIR           dune/grid           HINTS ${DUNE_ROOT}/include ${DUNE_ROOT}/dune-grid)
find_path(DUNE_LOCALFUNCTIONS_INCLUDE_DIR dune/localfunctions HINTS ${DUNE_ROOT}/include ${DUNE_ROOT}/dune-localfunctions)
find_path(DUNE_PDELAB_INCLUDE_DIR         dune/pdelab         HINTS ${DUNE_ROOT}/include ${DUNE_ROOT}/dune-pdelab)
find_path(DUNE_ISTL_INCLUDE_DIR           dune/istl           HINTS ${DUNE_ROOT}/include ${DUNE_ROOT}/dune-istl)
# find_path(DUNE_MULTIDOMAIN_INCLUDE_DIR    dune/grid/multidomaingrid           HINTS ${DUNE_ROOT}/include ${DUNE_ROOT}/dune-multidomaingrid)
# find_path(DUNE_TYPETREE_INCLUDE_DIR       dune/typetree       HINTS ${DUNE_ROOT}/include ${DUNE_ROOT}/dune-typetree)

find_library(DUNE_COMMON_LIBRARY         dunecommon   HINTS ${DUNE_ROOT}/lib ${DUNE_ROOT}/dune-common/lib/.libs)
find_library(DUNE_GRID_LIBRARY           dunegrid     HINTS ${DUNE_ROOT}/lib ${DUNE_ROOT}/dune-grid/lib/.libs)
find_library(DUNE_GEOMETRY_LIBRARY       dunegeometry HINTS ${DUNE_ROOT}/lib ${DUNE_ROOT}/dune-geometry/lib/.libs)
#find_library(DUNE_PDELAB_LIBRARY         dunepdelab   HINTS ${DUNE_ROOT}/lib ${DUNE_ROOT}/dune-pdelab/lib/.libs)

# there are no dune-istl and dune-localfunction libs
find_library(DUNE_LOCALFUNCTIONS_LIBRARY dunecommon   HINTS ${DUNE_ROOT}/lib ${DUNE_ROOT}/dune-pdelab/lib/.libs)
find_library(DUNE_ISTL_LIBRARY           dunecommon   HINTS ${DUNE_ROOT}/lib ${DUNE_ROOT}/dune-pdelab/lib/.libs)


set(DUNE_LIBRARIES    ${DUNE_GRID_LIBRARY} ${DUNE_GEOMETRY_LIBRARY} ${DUNE_COMMON_LIBRARY}) #${DUNE_PDELAB_LIBRARY} 
set(DUNE_INCLUDE_DIRS ${DUNE_COMMON_INCLUDE_DIR} ${DUNE_GEOMETRY_INCLUDE_DIR} ${DUNE_GRID_INCLUDE_DIR} ${DUNE_LOCALFUNCTIONS_INCLUDE_DIR} ${DUNE_ISTL_INCLUDE_DIR}) #${DUNE_PDELAB_INCLUDE_DIR} 

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(DUNE DEFAULT_MESSAGE DUNE_LIBRARIES DUNE_INCLUDE_DIRS)
find_package_handle_standard_args(Dune "dune-pdelab"         DUNE_PDELAB_LIBRARY 		DUNE_PDELAB_INCLUDE_DIR)
find_package_handle_standard_args(Dune "dune-grid"           DUNE_GRID_LIBRARY 		DUNE_GRID_INCLUDE_DIR)
find_package_handle_standard_args(Dune "dune-geometry"       DUNE_GEOMETRY_LIBRARY 		DUNE_GEOMETRY_INCLUDE_DIR)
find_package_handle_standard_args(Dune "dune-localfunctions" True		 	DUNE_LOCALFUNCTIONS_INCLUDE_DIR)
find_package_handle_standard_args(Dune "dune-istl"           TRUE           	DUNE_ISTL_INCLUDE_DIR)

