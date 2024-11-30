project(MatrixLib)

set(HEADERS
    include/Matrix.hpp
    include/SparseMatrix.hpp
    include/RowSparseMatrix.hpp
    include/PlotMatrix.hpp
)

set(SOURCES
    src/Matrix.cpp
    src/RowSparseMatrix.cpp
    src/PlotMatrix.cpp
)

add_library(MatrixLib SHARED ${SOURCES} ${HEADERS})

if (MPI_FOUND)
  include_directories(SYSTEM ${MPI_INCLUDE_PATH})
  target_link_libraries(MatrixLib ${MPI_LIBRARIES})
endif()

target_include_directories(MatrixLib PUBLIC include)