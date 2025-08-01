# Multivariate Newton-Raphson Solver in C

This project provides a C implementation of the multivariate Newton-Raphson method for finding the minimum of a function. It utilizes BLAS and LAPACK for efficient and robust linear algebra operations.

The solver is demonstrated with two examples within the `main` function of `newton_solver.c`:
1.  An optimization problem to fit a physical model of light intensity.
2.  A standard optimization test using the N-dimensional Rosenbrock function.

*Note: The current `main` function exits after running the light intensity model simulation and does not proceed to the Rosenbrock optimization.*

## Dependencies

To build and run this project, you will need:
- A C compiler (e.g., GCC, Clang)
- CMake (version 3.10 or later)
- BLAS and LAPACK libraries (e.g., OpenBLAS, MKL, or standard system libraries)

## Build Instructions

1.  **Clone the repository:**
    ```bash
    git clone <repository-url>
    cd <repository-directory>
    ```

2.  **Create a build directory:**
    ```bash
    mkdir build
    cd build
    ```

3.  **Run CMake to configure the project:**
    ```bash
    cmake ..
    ```
    This command will find the necessary libraries and generate the build files.

4.  **Compile the project:**
    ```bash
    make
    ```
    This will create an executable named `newton_solver` in the `build` directory.

## Usage

After building the project, you can run the executable from the project's root directory:
```bash
./build/newton_solver
```

The program will print the parameters of a randomly generated "true" intensity function, simulate some measurements, and then compute the initial distance for a test function.

## Code Structure

- **`newton_solver.c`**: The main source file containing:
    - Data structures for vectors (`Vector`) and matrices (`Matrix`).
    - A representation of the intensity function model (`Ifunc`, `Isample`).
    - The `newton_raphson` function, which performs the optimization.
    - The `user_function`, which defines the function to be minimized (currently the Rosenbrock function).
    - A `main` function that sets up and runs the intensity model fitting problem.
    - Helper functions for memory management and for solving linear systems using LAPACK's `dgesv` routine.

- **`CMakeLists.txt`**: The build script for CMake, which handles project configuration and links against BLAS and LAPACK.

- **`LICENSE`**: The MIT License file.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

Copyright (c) 2025 Olivier Guyon
