find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/contrib/include
                     ${CMAKE_CURRENT_SOURCE_DIR}/dune )

add_executable(test     test.cpp    )

target_link_libraries( test ${DUNE_LIBS} vtkHybrid vtkCommon vtkIO vtkFiltering )
