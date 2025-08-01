#ifndef NEWTON_SOLVER_H
#define NEWTON_SOLVER_H

#include "vector_matrix.h"

// Define a function pointer type for the objective function.
// This function should take a vector x and compute the value f(x),
// the gradient grad(f(x)), and the Hessian H(f(x)).
typedef double (*objective_function)(const Vector* x, Vector* grad, Matrix* hess);

void newton_raphson(Vector* x, objective_function func, int max_iter, double tolerance);

#endif // NEWTON_SOLVER_H
