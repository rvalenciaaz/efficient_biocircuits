#include <stdio.h>
#include <stdlib.h>
#include <cvode/cvode.h>             // prototypes for CVODE functions and constants
#include <nvector/nvector_serial.h>  // serial N_Vector types, functions, and macros
#include <sundials/sundials_math.h>  // definition of SUNRabs
#include <sunmatrix/sunmatrix_dense.h> // access to dense SUNMatrix
#include <sunlinsol/sunlinsol_dense.h> // access to dense SUNLinearSolver
#include <cvode/cvode_direct.h> // access to CVDls interface
#include <sundials/sundials_types.h>  // defs. of realtype, sunindextype

// Parameters for the simple gene expression model
#define BETA 1.0
#define GAMMA 0.5

// Function to compute the derivative
int simple_gene_expression(realtype t, N_Vector x, N_Vector xdot, void *user_data) {
    NV_Ith_S(xdot, 0) = BETA - GAMMA * NV_Ith_S(x, 0);
    return 0;
}

int main() {
    // Initial conditions
    realtype t0 = 0.0;
    realtype x0 = 0.0;
    realtype t;
    realtype dt = 0.1;
    realtype T = 50.0;

    // Create a serial vector of length 1 for the solution
    N_Vector x = N_VNew_Serial(1);
    NV_Ith_S(x, 0) = x0;

    // Create the CVODE memory block
    void *cvode_mem = CVodeCreate(CV_ADAMS, CV_NEWTON);
    if (cvode_mem == NULL) {
        fprintf(stderr, "Error in CVodeCreate\n");
        return 1;
    }

    // Initialize CVODE
    int flag = CVodeInit(cvode_mem, simple_gene_expression, t0, x);
    if (flag != CV_SUCCESS) {
        fprintf(stderr, "Error in CVodeInit\n");
        return 1;
    }

    // Specify the relative and absolute tolerances
    flag = CVodeSStolerances(cvode_mem, 1e-4, 1e-8);
    if (flag != CV_SUCCESS) {
        fprintf(stderr, "Error in CVodeSStolerances\n");
        return 1;
    }

    // Create the dense SUNMatrix
    SUNMatrix A = SUNDenseMatrix(1, 1);
    if (A == NULL) {
        fprintf(stderr, "Error in SUNDenseMatrix\n");
        return 1;
    }

    // Create the dense SUNLinearSolver
    SUNLinearSolver LS = SUNDenseLinearSolver(x, A);
    if (LS == NULL) {
        fprintf(stderr, "Error in SUNLinSol_Dense\n");
        return 1;
    }

    // Attach the linear solver to CVODE
    flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
    if (flag != CV_SUCCESS) {
        fprintf(stderr, "Error in CVodeSetLinearSolver\n");
        return 1;
    }

    // Open a file to save the results
    FILE *fp = fopen("simple_gene_expression_sundials.csv", "w");
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
            fprintf(stderr, "Error in CVode\n");
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
