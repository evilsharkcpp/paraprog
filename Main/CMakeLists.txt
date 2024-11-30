project(Main)

add_executable(Main main.cpp)

if(MPI_FOUND)
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})
    target_link_libraries(Main MatrixLib ${MPI_LIBRARIES})
else()
    target_link_libraries(Main MatrixLib)
endif()

