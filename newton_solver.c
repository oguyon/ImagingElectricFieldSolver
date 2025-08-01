#include "newton_solver.h"
#include "linear_solver.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/**
 * @brief Performs multivariate optimization using the Newton-Raphson method.
 */
void newton_raphson(Vector* x, objective_function obj_func, int max_iter, double tolerance) {
    int n = x->size;
    Vector* grad = create_vector(n);
    Matrix* hess = create_matrix(n, n);
    Vector* step = create_vector(n);

    for (int k = 0; k < max_iter; k++) {
        double f_val = obj_func(x, grad, hess);
        double grad_norm = 0.0;
        for (int i = 0; i < n; i++) {
            grad_norm += grad->data[i] * grad->data[i];
        }
        grad_norm = sqrt(grad_norm);

        printf("Iter %d: f(x)=%.6e, |grad|=%.6e\n", k, f_val, grad_norm);

        if (grad_norm < tolerance) {
            printf("Converged successfully!\n");
            break;
        }

        // Negate gradient for the linear system
        for (int i = 0; i < n; i++) {
            grad->data[i] = -grad->data[i];
        }

        if (solve_linear_system(hess, grad, step) != 0) {
            printf("Error: LAPACK solver failed. Matrix may be singular.\n");
            break;
        }

        // Update x
        for (int i = 0; i < n; i++) {
            x->data[i] += step->data[i];
        }

        if (k == max_iter - 1) {
            printf("Warning: Maximum iterations reached.\n");
        }
    }
    free_vector(grad);
    free_matrix(hess);
    free_vector(step);
}
