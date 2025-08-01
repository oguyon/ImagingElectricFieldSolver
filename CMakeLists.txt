# Specify the minimum version of CMake required to build this project.
cmake_minimum_required(VERSION 3.10)

# Define the project name.
project(NewtonSolverWithBLAS C)

# Set the C standard to C99.
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Find the BLAS and LAPACK libraries.
# CMake will find an available implementation, such as OpenBLAS, MKL, or the
# standard system libraries. The "REQUIRED" keyword will cause CMake to
# fail if it cannot find them.
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Add the executable target.
add_executable(newton_solver newton_solver.c)

# Link the found BLAS and LAPACK libraries to the executable.
# The ${LAPACK_LIBRARIES} variable holds the full set of libraries
# needed, including BLAS dependencies.
target_link_libraries(newton_solver PRIVATE ${LAPACK_LIBRARIES})

# Optional: Add an install rule.
install(TARGETS newton_solver DESTINATION bin)
