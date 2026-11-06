#include <stdio.h>
#include <math.h>
#include "matrix_print.h"

int QR_solver (int n, int m, double *A, double *x, double *b, int *p, double tol) {

    double modvv, modvb;

    #define A(i,j) A[i*n + j]

    // Checking if it's determined

    for (int i = 0; i < n; i++) if (fabs(A(i,p[i])) < tol) return -1;

    for (int k = 0; k < n; k++) {
        modvv = 1;
        modvb = b[k];
        for (int i = m-1; i > k; i--) {
            modvv += A(i,p[k]) * A(i,p[k]);
            modvb += A(i,p[k]) * b[i];
        }
        b[k] -= 2*modvb/modvv;
        for (int i = k+1; i < m; i++) {
            b[i] -= 2*modvb/modvv*A(i,p[k]);
        }
    }
    
    // Checking if it's compatible

    for (int i = n; i < m; i++) if (b[i] != 0) return -1;

    for (int i = n-1; i >= 0; i--) {
        x[p[i]] = b[i];
        for (int j = n-1; j > i; j--) {
            x[p[i]] -= A(i,p[j]) * x[p[j]];
        }
        x[p[i]] /= A(i,p[i]);
    }

    return 0;

}