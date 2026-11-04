#include <stdio.h>

int LU_solver(int n, double *A, int *p, double *x, double tol) {

    #define A(i,j) A[i*(n+1) + j]

    // We assume that A is of the LU form.


    double y[n];

    // Solving LUx=Pb. First, we solve Ly = Pb for y.

    for (int i = 0; i < n; i++) {
        y[i] = A(p[i],n);
        for (int j = 0; j < i; j++) y[i] -= A(p[i],j)*y[j];
    }

    // Now we solve Ux=y

    for (int i = n-1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = n-1; j > i; j--) x[i] -= A(p[i],j)*x[j];
        x[i] /= A(p[i],i);
    }

    return 0;
}