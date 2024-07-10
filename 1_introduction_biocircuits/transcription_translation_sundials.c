#include <stdio.h>
#include <stdlib.h>
#include <cvode/cvode.h>             // prototypes for CVODE functions and constants
#include <nvector/nvector_serial.h>  // serial N_Vector types, functions, and macros
#include <sundials/sundials_math.h>  // definition of SUNRabs
#include <sunmatrix/sunmatrix_dense.h> // access to dense SUNMatrix
#include <sunlinsol/sunlinsol_dense.h> // access to dense SUNLinearSolver
#include <cvode/cvode_direct.h> // access to CVDls interface
#include <sundials/sundials_types.h>  // defs. of realtype, sunindextype

// Parameters for the transcription and translation model
#define BETA_M 1.0
#define GAMMA_M 0.5
#define BETA_P 1.0
#define GAMMA_P 0.5

// Function to compute the derivatives
int transcription_translation(realtype t, N_Vector y, N_Vector ydot, void *user_data) {
    realtype m = NV_Ith_S(y, 0);
    realtype p = NV_Ith_S(y, 1);

    NV_Ith_S(ydot, 0) = BETA_M - GAMMA_M * m;
    NV_Ith_S(ydot, 1) = BETA_P * m - GAMMA_P * p;

    return 0;
}

int main() {
    // Initial conditions
    realtype t0 = 0.0;
    realtype m0 = 0.0;
    realtype p0 = 0.0;
    realtype t;
    realtype dt = 0.1;
    realtype T = 50.0;

    // Create a serial vector of length 2 for the solution
    N_Vector y = N_VNew_Serial(2);
    NV_Ith_S(y, 0) = m0;
    NV_Ith_S(y, 1) = p0;

    // Create the CVODE memory block
    void *cvode_mem = CVodeCreate(CV_ADAMS, CV_NEWTON);
    if (cvode_mem == NULL) {
        fprintf(stderr, "Error in CVodeCreate\n");
        return 1;
    }

    // Initialize CVODE
    int flag = CVodeInit(cvode_mem, transcription_translation, t0, y);
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
    SUNMatrix A = SUNDenseMatrix(2, 2);
    if (A == NULL) {
        fprintf(stderr, "Error in SUNDenseMatrix\n");
        return 1;
    }

    // Create the dense SUNLinearSolver
    SUNLinearSolver LS = SUNDenseLinearSolver(y, A);
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
    FILE *fp = fopen("transcription_translation_sundials.csv", "w");
    if (fp == NULL) {
        fprintf(stderr, "Error opening file!\n");
        return 1;
    }
    fprintf(fp, "Time,mRNA_concentration,Protein_concentration\n");

    // Time-stepping loop
    t = t0;
    while (t < T) {
        flag = CVode(cvode_mem, t + dt, y, &t, CV_NORMAL);
        if (flag != CV_SUCCESS) {
            fprintf(stderr, "Error in CVode\n");
            return 1;
        }
        fprintf(fp, "%f,%f,%f\n", t, NV_Ith_S(y, 0), NV_Ith_S(y, 1));
    }

    fclose(fp);

    // Free memory
    N_VDestroy(y);
    CVodeFree(&cvode_mem);
    SUNLinSolFree(LS);
    SUNMatDestroy(A);

    return 0;
}