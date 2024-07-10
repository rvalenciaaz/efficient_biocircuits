#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cvode/cvode.h>             // prototypes for CVODE functions and constants
#include <nvector/nvector_serial.h>  // serial N_Vector types, functions, and macros
#include <sunmatrix/sunmatrix_dense.h> // access to dense SUNMatrix
#include <sunlinsol/sunlinsol_dense.h> // access to dense SUNLinearSolver
#include <cvode/cvode_direct.h>      // access to CVDls interface
#include <sundials/sundials_types.h>  // definitions of realtype, sunindextype

// Model Parameters
#define PRODUCTION_RATE_X 0.1  // Example value for production rate of X
#define DEGRADATION_RATE_X 0.05 // Example value for degradation rate of X
#define PRODUCTION_RATE_Y 0.1  // Example value for production rate of Y
#define DEGRADATION_RATE_Y 0.05 // Example value for degradation rate of Y
#define HILL_COEFFICIENT 2     // Hill coefficient for non-linearity in response
#define COPY_NUMBER 1          // Baseline gene copy number

// Function to compute the derivatives of the IFFL system
int iffl_system(realtype t, N_Vector y, N_Vector ydot, void *user_data) {
    realtype x = NV_Ith_S(y, 0); // Concentration of X
    realtype y_conc = NV_Ith_S(y, 1); // Concentration of Y

    realtype dxdt = PRODUCTION_RATE_X * COPY_NUMBER - DEGRADATION_RATE_X * x;
    realtype dydt = (PRODUCTION_RATE_Y * COPY_NUMBER / (1 + pow(x, HILL_COEFFICIENT))) - DEGRADATION_RATE_Y * y_conc;

    NV_Ith_S(ydot, 0) = dxdt;
    NV_Ith_S(ydot, 1) = dydt;

    return 0;
}

// Main function to setup and solve the ODE
int main() {
    realtype t0 = 0.0, t = t0, T = 50.0, dt = 0.1;
    N_Vector y = N_VNew_Serial(2); // Vector for storing the concentrations of X and Y
    NV_Ith_S(y, 0) = 0.0; // Initial concentration of X
    NV_Ith_S(y, 1) = 0.0; // Initial concentration of Y

    // Create CVODE memory block and initialize solver
    void *cvode_mem = CVodeCreate(CV_ADAMS);
    CVodeInit(cvode_mem, iffl_system, t0, y);

    // Set scalar relative and absolute tolerances
    CVodeSStolerances(cvode_mem, 1e-4, 1e-8);

    // Create dense SUNMatrix for use in linear solvers
    SUNMatrix A = SUNDenseMatrix(2, 2);
    SUNLinearSolver LS = SUNDenseLinearSolver(y, A);
    CVDlsSetLinearSolver(cvode_mem, LS, A);

    // Open a file to save the results
    FILE *fp = fopen("iffl_simulation_results.csv", "w");
    if (fp == NULL) {
        fprintf(stderr, "Failed to open output file.\n");
        return 1;
    }
    fprintf(fp, "Time,X_Concentration,Y_Concentration\n");

    // Time-stepping loop
    while (t < T) {
        int flag = CVode(cvode_mem, t + dt, y, &t, CV_NORMAL);
        if (flag != CV_SUCCESS) {
            fprintf(stderr, "Solver error: %d\n", flag);
            break;
        }
        fprintf(fp, "%f,%f,%f\n", t, NV_Ith_S(y, 0), NV_Ith_S(y, 1));
    }

    // Close the file and free resources
    fclose(fp);
    N_VDestroy(y);
    CVodeFree(&cvode_mem);
    SUNLinSolFree(LS);
    SUNMatDestroy(A);

    return 0;
}
