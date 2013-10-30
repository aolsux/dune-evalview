#find_package(GooglePerfTools REQUIRED)
# find_package(VTK REQUIRED)
find_package(ALUGrid REQUIRED)
find_package(DUNE REQUIRED)
find_package(Metis REQUIRED)
find_package(MPI REQUIRED)
find_package(Boost COMPONENTS system filesystem serialization REQUIRED)

#include(${VTK_USE_FILE})


include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../contrib/include
                     ${CMAKE_CURRENT_SOURCE_DIR}/ )


include_directories(${Boost_INCLUDE_DIRS})
include_directories(${ALUGRID_INCLUDE_DIRS})
include_directories(${ALUGRID_INCLUDE_DIRS}/serial)
include_directories(${DUNE_INCLUDE_DIRS})
#include_directories(${METIS_INCLUDE_DIRS})
include_directories(${MPI_INCLUDE_PATH})
#include_directories(${GOOGLE_PERFTOOLS_INCLUDE_DIR})


add_executable(locate     locate.cpp    )

target_link_libraries( locate ${DUNE_LIBRARIES}
                            ${ALUGRID_LIBRARIES}
#                             ${METIS_LIBRARIES}
                            ${Boost_FILESYSTEM_LIBRARY}
                            ${Boost_SYSTEM_LIBRARY}
                            ${Boost_SERIALIZATION_LIBRARY}
                            ${MPI_LIBRARIES}
#                             ${VTK_LIBRARIES} 
#                             ${GOOGLE_PERFTOOLS_LIBRARIES}
                     )
                     
                     
add_executable(func     func.cpp    )

target_link_libraries( func ${DUNE_LIBRARIES}
                            ${ALUGRID_LIBRARIES}
                            ${METIS_LIBRARIES}
                            ${Boost_FILESYSTEM_LIBRARY}
                            ${Boost_SYSTEM_LIBRARY}
                            ${Boost_SERIALIZATION_LIBRARY}
                            ${MPI_LIBRARIES}
#                             ${VTK_LIBRARIES} 
#                             ${GOOGLE_PERFTOOLS_LIBRARIES}
                     )

