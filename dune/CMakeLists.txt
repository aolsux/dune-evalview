find_package(VTK REQUIRED)
find_package(ALUGrid REQUIRED)
find_package(DUNE REQUIRED)
find_package(Boost COMPONENTS system filesystem serialization REQUIRED)

include(${VTK_USE_FILE})


include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/contrib/include
                     ${CMAKE_CURRENT_SOURCE_DIR}/dune )


include_directories(${Boost_INCLUDE_DIRS})
include_directories(${DUNE_INCLUDE_DIRS})
include_directories(${ALUGRID_INCLUDE_DIRS})

add_executable(test     test.cpp    )

target_link_libraries( test ${DUNE_LIBRARIES}
                            ${ALUGRID_LIBRARIES}
                            ${METIS_LIBRARIES}
                            ${Boost_FILESYSTEM_LIBRARY}
                            ${Boost_SYSTEM_LIBRARY}
                            ${Boost_SERIALIZATION_LIBRARY}
                            vtkHybrid vtkCommon vtkIO vtkFiltering )

