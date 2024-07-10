#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define NUM_REACTIONS 8

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
double kap(double I) {
    return KAP_MAX * I / (I + KDA);
}

// Function for kout([RRp])
double kout(double RRp) {
    return KOUT_MAX * pow(RRp / KDR, N) / (pow(RRp / KDR, N) + 1);
}

// Propensity functions
void compute_propensities(double *a, double *x, double I) {
    double HK = x[0];
    double HKp = x[1];
    double RR = x[2];
    double RRp = x[3];
    double SR = x[4];
    double SRp = x[5];
    double PH = x[6];
    double Output = x[7];

    a[0] = BETA_HK; // Production of HK
    a[1] = DELTA * HK; // Degradation of HK
    a[2] = kap(I) * HK; // Autophosphorylation of HK
    a[3] = KT * HKp * RR; // Phosphotransfer from HKp to RR
    a[4] = KTC * HKp * SR; // Phosphotransfer from HKp to SR
    a[5] = kout(RRp); // Production of Output
    a[6] = DELTA * Output; // Degradation of Output
    a[7] = DELTA * HKp; // Degradation of HKp
}

// Perform one step of the Gillespie algorithm
void gillespie_step(double *x, double I, double *t) {
    double a[NUM_REACTIONS];
    compute_propensities(a, x, I);

    double a0 = 0.0;
    for (int i = 0; i < NUM_REACTIONS; i++) {
        a0 += a[i];
    }

    if (a0 == 0.0) return; // No more reactions

    // Determine the time until the next reaction
    double r1 = (double)rand() / RAND_MAX;
    double tau = -log(r1) / a0;
    *t += tau;

    // Determine which reaction occurs
    double r2 = (double)rand() / RAND_MAX;
    double sum = 0.0;
    int reaction = -1;
    for (int i = 0; i < NUM_REACTIONS; i++) {
        sum += a[i];
        if (r2 * a0 < sum) {
            reaction = i;
            break;
        }
    }

    // Update the state based on the selected reaction
    switch (reaction) {
        case 0: x[0] += 1; break; // Production of HK
        case 1: x[0] -= 1; break; // Degradation of HK
        case 2: x[0] -= 1; x[1] += 1; break; // Autophosphorylation of HK
        case 3: x[1] -= 1; x[3] += 1; break; // Phosphotransfer from HKp to RR
        case 4: x[1] -= 1; x[5] += 1; break; // Phosphotransfer from HKp to SR
        case 5: x[7] += 1; break; // Production of Output
        case 6: x[7] -= 1; break; // Degradation of Output
        case 7: x[1] -= 1; break; // Degradation of HKp
    }
}

// Function to solve the ODE using Gillespie algorithm and store results in an array
void solve_dichotomous_feedback(double *results, int n_steps, double dt, double I) {
    double t = 0.0;
    double x[8] = {0.0}; // Initial conditions: all concentrations start at 0

    for (int i = 0; i < n_steps; i++) {
        while (t < i * dt) {
            gillespie_step(x, I, &t);
        }
        for (int j = 0; j < 8; j++) {
            results[i * 8 + j] = x[j];
        }
    }
}
