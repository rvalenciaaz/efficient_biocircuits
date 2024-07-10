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

// Model parameters for gene expressions in a simple FFL
#define KXY 0.5
#define KXZ 0.5
#define KYZ 0.5
#define BETAY 1.0
#define BETAZ 1.0
#define GAMMAY 0.1
#define GAMMAZ 0.1
#define NXY 2
#define NXZ 2
#define NYZ 2

// Function to compute the derivatives for the FFL
int f(realtype t, N_Vector y, N_Vector ydot, void *user_data) {
    realtype Y = NV_Ith_S(y, 0);
    realtype Z = NV_Ith_S(y, 1);
    realtype X = *(realtype *)user_data; // X is provided as external input

    // Hill functions for activation
    realtype activation_XY = pow(X / KXY, NXY) / (1 + pow(X / KXY, NXY));
    realtype activation_XZ = pow(X / KXZ, NXZ) / (1 + pow(X / KXZ, NXZ));
    realtype activation_YZ = pow(Y / KYZ, NYZ) / (1 + pow(Y / KYZ, NYZ));

    // ODEs
    NV_Ith_S(ydot, 0) = BETAY * activation_XY - GAMMAY * Y;
    NV_Ith_S(ydot, 1) = BETAZ * (activation_XZ + activation_YZ) - GAMMAZ * Z;

    return 0;
}

int main() {
    realtype T = 10.0, t = 0.0, dt = 0.1;
    realtype X = 1.0; // Example constant input for X

    // Create a serial vector for storing Y and Z
    N_Vector y = N_VNew_Serial(2);
    NV_Ith_S(y, 0) = 0.0; // Initial Y
    NV_Ith_S(y, 1) = 0.0; // Initial Z

    // Create the CVODE memory block
    void *cvode_mem = CVodeCreate(CV_ADAMS, CV_NEWTON);
    if (cvode_mem == NULL) {
        fprintf(stderr, "Error in CVodeCreate\n");
        return 1;
    }

    // Initialize CVODE
    int flag = CVodeInit(cvode_mem, f, t, y);
    if (flag != CV_SUCCESS) {
        fprintf(stderr, "Error in CVodeInit\n");
        return 1;
    }

    // Set the user data to pass X to the derivative function
    CVodeSetUserData(cvode_mem, &X);

    // Specify the relative and absolute tolerances
    flag = CVodeSStolerances(cvode_mem, 1e-4, 1e-8);
    if (flag != CV_SUCCESS) {
        fprintf(stderr, "Error in CVodeSStolerances\n");
        return 1;
    }

    // Integrate over time
    while (t < T) {
        realtype tout = t + dt;
        flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
        if (flag != CV_SUCCESS) {
            fprintf(stderr, "Error in CVode at time %g\n", t);
            return 1;
        }
        printf("At time %g, Y = %g, Z = %g\n", t, NV_Ith_S(y, 0), NV_Ith_S(y, 1));
    }

    // Free resources
    N_VDestroy(y);
    CVodeFree(&cvode_mem);

    return 0;
}
