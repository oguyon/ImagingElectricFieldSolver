#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// --- Data Structures ---
// Represents a vector
typedef struct {
    int size;
    double* data;
} Vector;

// Represents a matrix
typedef struct {
    int rows;
    int cols;
    double** data;
} Matrix;



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





// --- LAPACK Function Prototype ---
// dgesv solves the system A*X = B for a general N-by-N matrix A.
// We need to declare it here to call it from C. The trailing underscore
// is a common convention for Fortran routines called from C.
extern void dgesv_(int* n, int* nrhs, double* a, int* lda, int* ipiv,
                   double* b, int* ldb, int* info);


// --- Function Prototypes ---
double user_function(const Vector* x, Vector* grad, Matrix* hess);
void newton_raphson(Vector* x0, int max_iter, double tolerance);
int solve_linear_system(Matrix* A, const Vector* b, Vector* x);
Vector* create_vector(int size);
Matrix* create_matrix(int size);
void free_vector(Vector* v);
void free_matrix(Matrix* m);




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






// --- Main Function ---
int main() {
    // number of WF modes, or poke modes
    int n_WFmodes = 2;
    // number of coherence channels, 2 for complex number monochromatic light
    int n_coch = 2;


    // Function parameters to be optimized
    int n_dims = 1 + n_WFmodes*(1+n_coch);
    printf("Number of dimensions for the optimization problem: %d\n", n_dims);
    Vector* x = create_vector(n_dims);
    for (int i = 0; i < n_dims; i++) {
        x->data[i] = 0.0;
    }

    // Define true intensity function
    printf("Creating Ifunction\n");
    Ifunc *ifc_true = create_Ifunc(n_coch, n_WFmodes);

    printf("Initializing Ifunc\n");
    // Intensity offset
    ifc_true->Iref = rand1();

    for (int ch=0; ch<n_coch; ch++)
    {
        // Amplitude offset for each channel
        ifc_true->Aref[ch] = rand1();

        // Amplitude coefficients
        for (int m=0; m<n_WFmodes; m++)
        {
            ifc_true->acoeff[ch][m] = rand1();
        }
    }


    // Define test intensity function
    printf("Creating Ifunction\n");
    Ifunc *ifc = create_Ifunc(n_coch, n_WFmodes);

    printf("Initializing Ifunc\n");
    double dampl = 0.1;
    // Intensity offset
    ifc->Iref = ifc_true->Iref + dampl * (1.0 - 2.0*rand1());

    for (int ch=0; ch<n_coch; ch++)
    {
        // Amplitude offset for each channel
        ifc->Aref[ch] = ifc_true->Aref[ch] + dampl * (1.0 - 2.0*rand1());

        // Amplitude coefficients
        for (int m=0; m<n_WFmodes; m++)
        {
            ifc->acoeff[ch][m] = ifc_true->acoeff[ch][m] + dampl * (1.0 - 2.0*rand1());
        }
    }





    printf("Bulding measurements\n");
    // Create simulated sample of measurements
    int NBsample = 10;

    Isample* isamplearr = create_Idamplearray(NBsample, n_WFmodes);

    for(int spi=0; spi<NBsample; spi++) {
        for(int wfm=0; wfm<n_WFmodes; wfm++) {
            isamplearr[spi].x[wfm] = 1.0 - 2.0*rand1();
        }
        evalIfunc(ifc_true, &isamplearr[spi]);
        printf("%4d  -> %f\n", spi, ifc_true->eval_I);

        // copy to measured instansity
        isamplearr[spi].Im = ifc_true->eval_I;
    }
    double distval = eval_dist_func(isamplearr, NBsample, ifc); 
    free_Isamplearray(isamplearr, NBsample);
    printf("DISTANCE = %f\n", distval);

    




    exit(0);



    printf("\nStarting Newton-Raphson Optimization for %d dimensions...\n", n_dims);
    printf("Initial guess: ");
    for(int i = 0; i < n_dims; i++) printf("%.4f ", x->data[i]);
    printf("\n");

    newton_raphson(x, 100, 1e-7);

    printf("----------------------------------------\n");
    printf("Optimization finished.\n");
    printf("Found minimum at: ");
    for(int i = 0; i < n_dims; i++) printf("%.6f ", x->data[i]);
    printf("\n");

    free_vector(x);
    return 0;
}


// --- Function Definitions ---


/**
 * @brief EXAMPLE: The N-dimensional Rosenbrock function.
 */
double user_function(const Vector* x, Vector* grad, Matrix* hess) {
    int n = x->size;
    double f_val = 0.0;

    for (int i = 0; i < n; i++) {
        grad->data[i] = 0.0;
        for (int j = 0; j < n; j++) {
            hess->data[i][j] = 0.0;
        }
    }

    for (int i = 0; i < n - 1; i++) {
        double x_i = x->data[i];
        double x_i_plus_1 = x->data[i+1];
        f_val += 100.0 * pow(x_i_plus_1 - x_i * x_i, 2) + pow(1.0 - x_i, 2);
        grad->data[i] += -400.0 * (x_i_plus_1 - x_i * x_i) * x_i - 2.0 * (1.0 - x_i);
        grad->data[i+1] += 200.0 * (x_i_plus_1 - x_i * x_i);
        hess->data[i][i] += 1200.0 * x_i * x_i - 400.0 * x_i_plus_1 + 2.0;
        hess->data[i+1][i+1] += 200.0;
        hess->data[i][i+1] += -400.0 * x_i;
        hess->data[i+1][i] += -400.0 * x_i;
    }
    return f_val;
}





/**
 * @brief Performs multivariate optimization using the Newton-Raphson method.
 */
void newton_raphson(Vector* x, int max_iter, double tolerance) {
    int n = x->size;
    Vector* grad = create_vector(n);
    Matrix* hess = create_matrix(n);
    Vector* step = create_vector(n);

    for (int k = 0; k < max_iter; k++) {
        double f_val = user_function(x, grad, hess);
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

        for (int i = 0; i < n; i++) {
            grad->data[i] = -grad->data[i];
        }

        if (solve_linear_system(hess, grad, step) != 0) {
            printf("Error: LAPACK solver failed. Matrix may be singular.\n");
            break;
        }

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


// --- Memory Management Helpers ---
Vector* create_vector(int size) {
    Vector* v = (Vector*)malloc(sizeof(Vector));
    v->size = size;
    v->data = (double*)malloc(size * sizeof(double));
    return v;
}

Matrix* create_matrix(int size) {
    Matrix* m = (Matrix*)malloc(sizeof(Matrix));
    m->rows = size;
    m->cols = size;
    m->data = (double**)malloc(size * sizeof(double*));
    for (int i = 0; i < size; i++) {
        m->data[i] = (double*)malloc(size * sizeof(double));
    }
    return m;
}

void free_vector(Vector* v) {
    if (v) { free(v->data); free(v); }
}

void free_matrix(Matrix* m) {
    if (m) {
        for (int i = 0; i < m->rows; i++) { free(m->data[i]); }
        free(m->data);
        free(m);
    }
}
