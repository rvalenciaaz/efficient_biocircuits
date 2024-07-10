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

// Parameters for the model
#define BETA_HK 1.0
#define BETA_RR 1.0
#define BETA_SR 1.0
#define BETA_PH 1.0
#define DELTA 0.1
#define KAP_MAX 1.0
#define KDA 1.0
#define KT 0.1
#define KTC 0.1
#define KP 0.1
#define KPC 0.1
#define KOUT_MAX 1.0
#define KDR 1.0
#define N 2

// Function for kap(I)
realtype kap(realtype I) {
    return KAP_MAX * I / (I + KDA);
}

// Function for kout([RRp])
realtype kout(realtype RRp) {
    return KOUT_MAX * pow(RRp / KDR, N) / (pow(RRp / KDR, N) + 1);
}

// Function to compute the derivatives
int dichotomous_feedback(realtype t, N_Vector y, N_Vector ydot, void *user_data) {
    realtype *params = (realtype *)user_data;
    realtype I = params[0];

    realtype HK = NV_Ith_S(y, 0);
    realtype HKp = NV_Ith_S(y, 1);
    realtype RR = NV_Ith_S(y, 2);
    realtype RRp = NV_Ith_S(y, 3);
    realtype SR = NV_Ith_S(y, 4);
    realtype SRp = NV_Ith_S(y, 5);
    realtype PH = NV_Ith_S(y, 6);
    realtype Output = NV_Ith_S(y, 7);

    NV_Ith_S(ydot, 0) = BETA_HK - DELTA * HK - kap(I) * HK + KT * HKp * RR + KTC * HKp * SR;
    NV_Ith_S(ydot, 1) = -KT * HKp * RR + kap(I) * HK - DELTA * HKp - KTC * HKp * SR;
    NV_Ith_S(ydot, 2) = BETA_RR - DELTA * RR - KT * HKp * RR + KP * HK * RRp + KPC * PH * RRp;
    NV_Ith_S(ydot, 3) = -DELTA * RRp + KT * HKp * RR - KP * HK * RRp - KPC * PH * RRp;
    NV_Ith_S(ydot, 4) = BETA_SR - DELTA * SR - KTC * HKp * SR + KPC * HK * SRp;
    NV_Ith_S(ydot, 5) = -DELTA * SRp + KTC * HKp * SR - KPC * HK * SRp;
    NV_Ith_S(ydot, 6) = BETA_PH - DELTA * PH;
    NV_Ith_S(ydot, 7) = kout(RRp) - DELTA * Output;

    return 0;
}

// Function to solve the ODE and store results in an array
void solve_dichotomous_feedback(double *results, int n_steps, double dt, double I) {
    realtype t0 = 0.0;
    realtype t;
    realtype T = dt * (n_steps - 1);

    // Initial conditions
    N_Vector y = N_VNew_Serial(8);
    NV_Ith_S(y, 0) = 0.0;  // HK
    NV_Ith_S(y, 1) = 0.0;  // HKp
    NV_Ith_S(y, 2) = 0.0;  // RR
    NV_Ith_S(y, 3) = 0.0;  // RRp
    NV_Ith_S(y, 4) = 0.0;  // SR
    NV_Ith_S(y, 5) = 0.0;  // SRp
    NV_Ith_S(y, 6) = 0.0;  // PH
    NV_Ith_S(y, 7) = 0.0;  // Output

    // Parameters
    realtype params[1] = {I};

    // Create the CVODE memory block
    void *cvode_mem = CVodeCreate(CV_ADAMS, CV_NEWTON);
    if (cvode_mem == NULL) {
        fprintf(stderr, "Error in CVodeCreate\n");
        return;
    }

    // Initialize CVODE
    int flag = CVodeInit(cvode_mem, dichotomous_feedback, t0, y);
    if (flag != CV_SUCCESS) {
        fprintf(stderr, "Error in CVodeInit\n");
        return;
    }

    // Specify the relative and absolute tolerances
    flag = CVodeSStolerances(cvode_mem, 1e-4, 1e-8);
    if (flag != CV_SUCCESS) {
        fprintf(stderr, "Error in CVodeSStolerances\n");
        return;
    }

    // Create the dense SUNMatrix
    SUNMatrix A = SUNDenseMatrix(8, 8);
    if (A == NULL) {
        fprintf(stderr, "Error in SUNDenseMatrix\n");
        return;
    }

    // Create the dense SUNLinearSolver
    SUNLinearSolver LS = SUNDenseLinearSolver(y, A);
    if (LS == NULL) {
        fprintf(stderr, "Error in SUNLinSol_Dense\n");
        return;
    }

    // Attach the linear solver to CVODE
    flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
    if (flag != CV_SUCCESS) {
        fprintf(stderr, "Error in CVodeSetLinearSolver\n");
        return;
    }

    // Set the user data
    flag = CVodeSetUserData(cvode_mem, params);
    if (flag != CV_SUCCESS) {
        fprintf(stderr, "Error in CVodeSetUserData\n");
        return;
    }

    // Time-stepping loop
    t = t0;
    for (int i = 0; i < n_steps; i++) {
        flag = CVode(cvode_mem, t + dt, y, &t, CV_NORMAL);
        if (flag != CV_SUCCESS) {
            fprintf(stderr, "Error in CVode\n");
            return;
        }
        for (int j = 0; j < 8; j++) {
            results[i * 8 + j] = NV_Ith_S(y, j);
        }
    }

    // Free
    N_VDestroy(y);
    CVodeFree(&cvode_mem);
    SUNLinSolFree(LS);
    SUNMatDestroy(A);
}
