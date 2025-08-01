#ifndef IFUNC_H
#define IFUNC_H

#include "vector_matrix.h"

// --- Data Structures ---
// Parameters of function I (single sample intensity)
typedef struct {
    // number of coherence channels
    int n_coch;
    // number of WF modes
    int n_WFmode;

    // total number of function parameters
    int n_dims;

    // Intensity offset
    double Iref;
    // Amplitude offset (per channel)
    double* Aref;
    // Amplitude coefficients
    double** acoeff;

    // Computed I, its first and second order derivatives
    double eval_I;
    double* eval_dI;
    double** eval_ddI;
} Ifunc;


// Sample point
typedef struct {
    // Coordinates
    double *x;
    // Measured intensity
    double Im;
} Isample;

// --- Function Prototypes ---
Isample* create_Idamplearray(int n_sample, int n_WFmode);
void free_Isamplearray(Isample* isample, int n_sample);
Ifunc* create_Ifunc(int n_coch, int n_WFmode);
void free_Ifunc(Ifunc* ifc);
void Ifunc_2_optvars(Ifunc *ifc, double *optvars);
void optvars_2_Ifunc(Ifunc *ifc, double *optvars);
double evalIfunc(Ifunc* ifc, Isample* isample);
double eval_dist_func(Isample* isamplearr, int NBsample, Ifunc* ifc, double *regalpha, double *regvec);
double rand1();

#endif // IFUNC_H
