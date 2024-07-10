#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cvode/cvode.h>             // prototypes for CVODE functions and constants
#include <nvector/nvector_serial.h>  // serial N_Vector types, functions, and macros
#include <sundials/sundials_math.h>  // definition of SUNRabs
#include <sunmatrix/sunmatrix_dense.h> // access to dense SUNMatrix
#include <sunlinsol/sunlinsol_dense.h> // access to dense SUNLinearSolver
#include <cvode/cvode_direct.h> // access to CVDls interface
#include <sundials/sundials_types.h>  // defs. of realtype, sunindextype

// Parameters for the autoregulatory gene expression model
#define BETA 10.0  // Maximum production rate
#define GAMMA 1.0  // Degradation rate
#define K 3.0      // Hill constant
#define N 5        // Hill coefficient

// Function to compute the derivative
int autoregulatory_gene_expression(realtype t, N_Vector x, N_Vector xdot, void *user_data) {
    realtype x_val = NV_Ith_S(x, 0);  // Get the current value of x
    NV_Ith_S(xdot, 0) = BETA * pow(x_val, N) / (pow(K, N) + pow(x_val, N)) - GAMMA * x_val;
    return 0;
}

int main() {
    // Initial conditions
    realtype t0 = 0.0;
    realtype x0 = 0.1;  // Initial concentration of the protein
    realtype t;
    realtype dt = 0.1;
    realtype T = 50.0;  // Total time for simulation

    // Create a serial vector of length 1 for the solution
    N_Vector x = N_VNew_Serial(1);
    if (x == NULL) {
        fprintf(stderr, "Failed to create a new vector for x.\n");
        return 1;
    }
    NV_Ith_S(x, 0) = x0;

    // Create the CVODE memory block
    void *cvode_mem = CVodeCreate(CV_ADAMS);
    if (cvode_mem == NULL) {
        fprintf(stderr, "Error in CVodeCreate\n");
        return 1;
    }

    // Initialize CVODE with the user's function, initial time t0, and initial dependent variable vector x
    int flag = CVodeInit(cvode_mem, autoregulatory_gene_expression, t0, x);
    if (flag != CV_SUCCESS) {
        fprintf(stderr, "Error in CVodeInit: %d\n", flag);
        return 1;
    }

    // Set scalar relative and absolute tolerances
    flag = CVodeSStolerances(cvode_mem, 1e-4, 1e-8);
    if (flag != CV_SUCCESS) {
        fprintf(stderr, "Error in CVodeSStolerances: %d\n", flag);
        return 1;
    }

    // Create the dense SUNMatrix for use in linear solves
    SUNMatrix A = SUNDenseMatrix(1, 1);
    if (A == NULL) {
        fprintf(stderr, "Error creating dense SUNMatrix.\n");
        return 1;
    }

    // Create the dense SUNLinearSolver object
    SUNLinearSolver LS = SUNDenseLinearSolver(x, A);
    if (LS == NULL) {
        fprintf(stderr, "Error creating dense SUNLinearSolver.\n");
        return 1;
    }

    // Attach the linear solver
    flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
    if (flag != CV_SUCCESS) {
        fprintf(stderr, "Error in CVDlsSetLinearSolver: %d\n", flag);
        return 1;
    }

    // Open a file to save the results
    FILE *fp = fopen("autoregulatory_gene_expression.csv", "w");
    if (fp == NULL) {
        fprintf(stderr, "Error opening file!\n");
        return 1;
    }
    fprintf(fp, "Time,Protein_concentration\n");

    // Time-stepping loop
    t = t0;
    while (t < T) {
        flag = CVode(cvode_mem, t + dt, x, &t, CV_NORMAL);
        if (flag != CV_SUCCESS) {
            fprintf(stderr, "Error in CVode: %d\n", flag);
            return 1;
        }
        fprintf(fp, "%f,%f\n", t, NV_Ith_S(x, 0));
    }

    fclose(fp);

    // Free memory
    N_VDestroy(x);
    CVodeFree(&cvode_mem);
    SUNLinSolFree(LS);
    SUNMatDestroy(A);

    return 0;
}
