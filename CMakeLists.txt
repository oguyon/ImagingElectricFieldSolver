# Specify the minimum version of CMake required to build this project.
cmake_minimum_required(VERSION 3.10)

# Define the project name.
project(NewtonSolverWithBLAS C)

# Set the C standard to C99.
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Find the BLAS and LAPACK libraries.
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Add the executable target with all source files.
add_executable(newton_solver
    main.c
    vector_matrix.c
    linear_solver.c
    ifunc.c
    newton_solver.c
)

# Link the found BLAS and LAPACK libraries to the executable.
target_link_libraries(newton_solver PRIVATE ${LAPACK_LIBRARIES})

# Optional: Add an install rule.
install(TARGETS newton_solver DESTINATION bin)
