#include <stdio.h>

// Parameters for the simple gene expression model
#define BETA 1.0
#define GAMMA 0.5

// Function to compute the derivative
double simple_gene_expression(double x) {
    return BETA - GAMMA * x;
}

int main() {
    double x = 0; // initial condition
    double t = 0;
    double dt = 0.1; // time step
    double T = 50; // total time

    FILE *fp = fopen("simple_gene_expression.csv", "w");
    if (fp == NULL) {
        printf("Error opening file!\n");
        return 1;
    }

    fprintf(fp, "Time,Protein_concentration\n");

    while (t <= T) {
        fprintf(fp, "%f,%f\n", t, x);
        x += dt * simple_gene_expression(x);
        t += dt;
    }

    fclose(fp);
    return 0;
}
