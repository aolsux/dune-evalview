project(dune-evalview)

cmake_minimum_required(VERSION 2.8)

set(CMAKE_BUILD_TYPE "Release" )

set(DEFINES "-Wno-unknown-pragmas -Wdelete-non-virtual-dtor -Wno-deprecated-declarations")

set(CMAKE_C_FLAGS_RELWITHDEBINFO "${DEFINES} -march=native -fopenmp -std=c++11 -g2 -O0 -Wall -Wno-long-long -pipe -fPIC -malign-double -pedantic")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${DEFINES} -march=native -fopenmp -std=c++11 -g2 -O0 -Wall -Wno-long-long -pipe -fPIC -malign-double -pedantic")

set(CMAKE_C_FLAGS_DEBUG "${DEFINES} -march=native -fopenmp -std=c++11 -g3 -O0 -Wall -Wno-long-long -pipe -fPIC -malign-double -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "${DEFINES} -march=native -fopenmp -std=c++11 -g3 -O0 -Wall -Wno-long-long -pipe -fPIC -malign-double -pedantic")

set(CMAKE_C_FLAGS_RELEASE "${DEFINES} -march=native -fopenmp -std=c++11 -O3 -g3 -Wall -Wno-long-long -pipe -fPIC -malign-double -pedantic" )
set(CMAKE_CXX_FLAGS_RELEASE "${DEFINES} -march=native -fopenmp -std=c++11 -O3 -g3 -Wall -Wno-long-long -pipe -fPIC -malign-double -pedantic")

set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--export-dynamic,--no-undefined,--as-needed")

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH} )

add_subdirectory( dune )
