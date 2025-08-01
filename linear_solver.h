#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H

#include "vector_matrix.h"

// Solves the linear system A*x = b using LAPACK
int solve_linear_system(Matrix* A, const Vector* b, Vector* x);

#endif // LINEAR_SOLVER_H
