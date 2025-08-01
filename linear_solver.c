#include "linear_solver.h"
#include <stdlib.h>

// --- LAPACK Function Prototype ---
// dgesv solves the system A*X = B for a general N-by-N matrix A.
// We need to declare it here to call it from C. The trailing underscore
// is a common convention for Fortran routines called from C.
extern void dgesv_(int* n, int* nrhs, double* a, int* lda, int* ipiv,
                   double* b, int* ldb, int* info);

/**
 * @brief Solves the linear system Ax = b using LAPACK's dgesv routine.
 */
int solve_linear_system(Matrix* A, const Vector* b, Vector* x) {
    int n = A->rows;
    int nrhs = 1; // Number of right-hand sides
    int lda = n;
    int ldb = n;
    int info;

    // LAPACK overwrites input arrays, so we must use copies.
    // LAPACK also expects column-major order, so we transpose A.
    double* a_col_major = (double*)malloc(n * n * sizeof(double));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            a_col_major[j * n + i] = A->data[i][j]; // Transpose
        }
    }

    // Copy b to x, as dgesv will overwrite b with the solution.
    for (int i = 0; i < n; i++) {
        x->data[i] = b->data[i];
    }

    // Allocate pivot array for dgesv
    int* ipiv = (int*)malloc(n * sizeof(int));

    // Call the LAPACK routine
    dgesv_(&n, &nrhs, a_col_major, &lda, ipiv, x->data, &ldb, &info);

    // Cleanup
    free(a_col_major);
    free(ipiv);

    // info=0 means success. info > 0 means matrix is singular.
    return info == 0 ? 0 : -1;
}
