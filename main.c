#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ifunc.h"
#include "newton_solver.h"
#include "vector_matrix.h"


// --- Global variables for the objective function ---
// This is not ideal, but for this specific problem structure,
// it's a simple way to pass extra data to the objective function.
static Isample* g_isamplearr = NULL;
static int g_NBsample = 0;
static Ifunc* g_ifc = NULL;

/**
 * @brief Objective function for the Ifunc optimization problem.
 *
 * This function calculates the sum of squared differences between measured
 * and modeled intensities. It serves as the objective function for the
 * Newton-Raphson solver.
 */
double ifunc_objective_function(const Vector* x, Vector* grad, Matrix* hess) {
    // Map the optimization variables (from vector x) to the Ifunc structure
    optvars_2_Ifunc(g_ifc, x->data);

    // Initialize gradient and hessian to zero
    for (int i = 0; i < x->size; i++) {
        grad->data[i] = 0.0;
        for (int j = 0; j < x->size; j++) {
            hess->data[i][j] = 0.0;
        }
    }

    double total_dist = 0.0;

    // Accumulate distance, gradient, and hessian over all samples
    for (int i = 0; i < g_NBsample; i++) {
        evalIfunc(g_ifc, &g_isamplearr[i]);
        double diff = g_isamplearr[i].Im - g_ifc->eval_I;
        total_dist += diff * diff;

        // Accumulate gradient: -2 * diff * dI/dp
        for (int j = 0; j < x->size; j++) {
            grad->data[j] += -2.0 * diff * g_ifc->eval_dI[j];
        }

        // Accumulate Hessian: 2 * (dI/dp_j * dI/dp_k - diff * d2I/dp_j dp_k)
        for (int j = 0; j < x->size; j++) {
            for (int k = 0; k < x->size; k++) {
                hess->data[j][k] += 2.0 * (g_ifc->eval_dI[j] * g_ifc->eval_dI[k] - diff * g_ifc->eval_ddI[j][k]);
            }
        }
    }

    return total_dist;
}


// --- Main Function ---
int main() {
    // number of WF modes, or poke modes
    int n_WFmodes = 2;
    // number of coherence channels, 2 for complex number monochromatic light
    int n_coch = 2;

    // Function parameters to be optimized
    int n_dims = 1 + n_WFmodes * (1 + n_coch);
    printf("Number of dimensions for the optimization problem: %d\n", n_dims);

    // --- Ground Truth Simulation ---
    printf("Creating true Ifunction\n");
    Ifunc *ifc_true = create_Ifunc(n_coch, n_WFmodes);

    printf("Initializing true Ifunc with random values\n");
    ifc_true->Iref = rand1();
    for (int ch = 0; ch < n_coch; ch++) {
        ifc_true->Aref[ch] = rand1();
        for (int m = 0; m < n_WFmodes; m++) {
            ifc_true->acoeff[ch][m] = rand1();
        }
    }

    // --- Measurement Simulation ---
    printf("Generating simulated measurements\n");
    int NBsample = 10;
    Isample* isamplearr = create_Idamplearray(NBsample, n_WFmodes);
    for (int spi = 0; spi < NBsample; spi++) {
        for (int wfm = 0; wfm < n_WFmodes; wfm++) {
            isamplearr[spi].x[wfm] = 1.0 - 2.0 * rand1();
        }
        isamplearr[spi].Im = evalIfunc(ifc_true, &isamplearr[spi]);
        printf("  Sample %d: Measured Intensity = %f\n", spi, isamplearr[spi].Im);
    }

    // --- Initial Guess ---
    printf("\nCreating initial guess Ifunc\n");
    Ifunc *ifc_guess = create_Ifunc(n_coch, n_WFmodes);
    double dampl = 0.1;
    ifc_guess->Iref = ifc_true->Iref + dampl * (1.0 - 2.0 * rand1());
    for (int ch = 0; ch < n_coch; ch++) {
        ifc_guess->Aref[ch] = ifc_true->Aref[ch] + dampl * (1.0 - 2.0 * rand1());
        for (int m = 0; m < n_WFmodes; m++) {
            ifc_guess->acoeff[ch][m] = ifc_true->acoeff[ch][m] + dampl * (1.0 - 2.0 * rand1());
        }
    }

    // --- Setup for Optimization ---
    // Set global variables to be used by the objective function
    g_isamplearr = isamplearr;
    g_NBsample = NBsample;
    g_ifc = ifc_guess;

    // Create the vector of optimization variables from the initial guess
    Vector* x = create_vector(n_dims);
    Ifunc_2_optvars(ifc_guess, x->data);

    printf("\nStarting Newton-Raphson Optimization for %d dimensions...\n", n_dims);
    printf("Initial guess: ");
    for(int i = 0; i < n_dims; i++) printf("%.4f ", x->data[i]);
    printf("\n");

    // --- Run Optimization ---
    newton_raphson(x, ifunc_objective_function, 100, 1e-7);

    printf("----------------------------------------\n");
    printf("Optimization finished.\n");

    // --- Results ---
    // Map final optimization variables back to Ifunc structure
    optvars_2_Ifunc(ifc_guess, x->data);
    printf("Found minimum parameters:\n");
    printf("  Iref: True=%.6f, Found=%.6f\n", ifc_true->Iref, ifc_guess->Iref);
    for(int ch=0; ch<n_coch; ch++) {
        printf("  Channel %d:\n", ch);
        printf("    Aref: True=%.6f, Found=%.6f\n", ifc_true->Aref[ch], ifc_guess->Aref[ch]);
        for (int m=0; m<n_WFmodes; m++) {
            printf("    acoeff[%d]: True=%.6f, Found=%.6f\n", m, ifc_true->acoeff[ch][m], ifc_guess->acoeff[ch][m]);
        }
    }

    // --- Cleanup ---
    free_vector(x);
    free_Ifunc(ifc_true);
    free_Ifunc(ifc_guess);
    free_Isamplearray(isamplearr, NBsample);

    return 0;
}
