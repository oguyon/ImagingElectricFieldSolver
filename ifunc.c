#include "ifunc.h"
#include <stdlib.h>
#include <stdio.h>

// Make a random number between 0.0 and 1.0
double rand1()
{
    return (double)rand() / (double)RAND_MAX;
}

Isample* create_Idamplearray(int n_sample, int n_WFmode)
{
    Isample* isample = (Isample*) malloc(sizeof(Isample)*n_sample);
    for(int si=0; si<n_sample; si++){
        isample[si].x = (double*) malloc(sizeof(double)*n_WFmode);
        isample[si].Im = 0.0;
    }
    return isample;
}

void free_Isamplearray(Isample* isample, int n_sample)
{
    for(int si=0; si<n_sample; si++){
        free(isample[si].x);
    }
    // There is a memory leak in the original code.
    // The array holding the samples also needs to be freed.
    free(isample);
}


Ifunc* create_Ifunc(int n_coch, int n_WFmode)
{
    Ifunc* ifc = (Ifunc*) malloc(sizeof(Ifunc));
    ifc->n_coch = n_coch;
    ifc->n_WFmode = n_WFmode;
    ifc->Aref = (double*) malloc(sizeof(double)*n_coch);
    ifc->acoeff = (double**) malloc(sizeof(double*)*n_coch);
    for (int cch = 0; cch < n_coch; cch++) {
        ifc->acoeff[cch]  = (double*)malloc(n_WFmode * sizeof(double));
    }

    ifc->n_dims = 1 + ifc->n_WFmode*(1 + ifc->n_coch);
    ifc->eval_dI = (double*) malloc(sizeof(double)*ifc->n_dims);
    ifc->eval_ddI = (double**) malloc(sizeof(double*)*ifc->n_dims);
    for (int dim = 0; dim < ifc->n_dims; dim++) {
        ifc->eval_ddI[dim]  = (double*)malloc(ifc->n_dims * sizeof(double));
    }

    return ifc;
}

void free_Ifunc(Ifunc* ifc)
{
    free(ifc->Aref);
    for (int cch = 0; cch < ifc->n_coch; cch++) {
        free(ifc->acoeff[cch]);
    }
    free(ifc->acoeff);

    int n_dims = 1 + ifc->n_WFmode*(1 + ifc->n_coch);
    free(ifc->eval_dI);
    for (int dim = 0; dim < n_dims; dim++) {
        free(ifc->eval_ddI[dim]);
    }
    free(ifc->eval_ddI);

    free(ifc);
}


// map I function parameters to optimization variables
void Ifunc_2_optvars( Ifunc *ifc, double *optvars)
{
    optvars[0] = ifc->Iref;
    for (int cch = 0; cch < ifc->n_coch; cch++) {
        optvars[1 + cch*(1+ifc->n_WFmode)] = ifc->Aref[cch];
        for (int m=0; m<ifc->n_WFmode; m++) {
            optvars[2+cch*(1+ifc->n_WFmode) + m] = ifc->acoeff[cch][m];
        }
    }
}

// map optimization variables to I function
void optvars_2_Ifunc( Ifunc *ifc, double *optvars)
{
    ifc->Iref = optvars[0];
    for (int cch = 0; cch < ifc->n_coch; cch++) {
        ifc->Aref[cch] = optvars[1 + cch*(1+ifc->n_WFmode)];
        for (int m=0; m<ifc->n_WFmode; m++) {
            ifc->acoeff[cch][m] = optvars[2+cch*(1+ifc->n_WFmode) + m];
        }
    }
}


// Evaluate intensity function, its 1st and 2nd order derivatives
// results stored in ifc
double evalIfunc(Ifunc* ifc, Isample* isample)
{
    double val = 0.0;

    // initialize 1st and 2nd derivatives to zero
    for(int i=0; i < ifc->n_dims; i++){
        ifc->eval_dI[i] = 0.0;
        for(int j=0; j < ifc->n_dims; j++){
            ifc->eval_ddI[i][j] = 0.0;
        }
    }

    // amplitude for each channel
    double* amp = (double*) malloc(sizeof(double)*ifc->n_coch);

    val += ifc->Iref;
    ifc->eval_dI[0] = 1.0;
    ifc->eval_ddI[0][0] = 0.0;
    for(int i=1; i < ifc->n_dims; i++){
        ifc->eval_ddI[i][0] = 0.0;
        ifc->eval_ddI[0][i] = 0.0;
    }

    for(int cch=0; cch<ifc->n_coch; cch++)
    {
        amp[cch] = ifc->Aref[cch];
        for(int m=0; m < ifc->n_WFmode; m++) {
            amp[cch] += isample->x[m] * ifc->acoeff[cch][m];
        }
        val += amp[cch]*amp[cch];

        // dI / dAref(ch)
        ifc->eval_dI[1 + (1+ifc->n_WFmode)*cch ] = 2.0 * amp[cch];

        // dI / da(m,cch)
         for(int m=0; m < ifc->n_WFmode; m++) {
            ifc->eval_dI[1 + (1+ifc->n_WFmode)*cch + 1+m] = 2.0 * isample->x[m] * amp[cch];
         }
    }

    // 2nd order derivatives

    // ddI / da(m0,cch0) / da(m1,cch1)
    // is zero if cch0 != cch1
    // so we just use cch index
    for(int cch=0; cch<ifc->n_coch; cch++) {
        for(int m0=0; m0<ifc->n_WFmode; m0++) {
            int i0 = 1 + (1 + ifc->n_WFmode)*cch + 1+m0;
            for(int m1=0; m1<ifc->n_WFmode; m1++) {
                int i1 = 1 + (1 + ifc->n_WFmode)*cch + 1+m1;
                ifc->eval_ddI[i0][i1] = 2.0 * isample->x[i0] * isample->x[i1];
            }
        }
    }

    // ddI / da(m0,cch0) / dAref(cch1)
    // ddI / dAref(cch0) / da(m1,cch1)
    // we recall that dI/da(m,cch0) = 2 xi amp[cch0]
    // and damp[cch0]/dAref[cch1] = 1 if cch0 = cch1, zero otherwise
    // so again we use a single cch index
    for(int cch=0; cch<ifc->n_coch; cch++) {
        for(int m=0; m<ifc->n_WFmode; m++) {
            // index i0 points to da(m,cch)
            int i0 = 1 + (1 + ifc->n_WFmode)*cch + 1+m;
            // index i1 points to dAref(cch)
            int i1 = 1 + (1 + ifc->n_WFmode)*cch;
            ifc->eval_ddI[i0][i1] = 2.0 * isample->x[m];
            ifc->eval_ddI[i1][i0] = 2.0 * isample->x[m];
        }
    }

    // ddI / dAref(cch0) / dAref(cch1)
    // we recall that dI/dAref(cch0) = 2.0 * amp[cch0]
    // and damp[cch0]/dAref[cch1] = 1 if cch0 = cch1, zero otherwise
    // so ddI = 2.0 if cch0 = cch1, 0 otherwise
    for(int cch=0; cch<ifc->n_coch; cch++) {
        int i = 1 + (1 + ifc->n_WFmode)*cch;
        ifc->eval_ddI[i][i] = 2.0;
    }


    ifc->eval_I = val;

    free(amp);

    return val;
}



// regalpha is contraining the solution to stay close to zero
double eval_dist_func(
    Isample* isamplearr, int NBsample,
    Ifunc* ifc,
    double *regalpha, double *regvec
)
{
    double distval = 0.0;
    for(int spi=0; spi<NBsample; spi++) {
        evalIfunc(ifc, &isamplearr[spi]);
        double Idiff = (isamplearr[spi].Im - ifc->eval_I);

        distval += Idiff*Idiff;
    }

    // Add regularization


    return distval;
}
