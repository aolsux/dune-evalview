find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(test     test.cpp    )

target_link_libraries( hang ${DUNE_LIBS} vtkHybrid vtkCommon vtkIO vtkFiltering )

