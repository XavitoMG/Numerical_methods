#include <math.h>
#include <stdio.h>
#include "matrix_print.h"

int gauss_seidel (int n, double *A, double *x, double *b, double tol, int stepmax) {

    #define A(i,j) A[i*n + j]

    double sum;
    double error = 0;
    double maxerror = 0;

    for (int i = 0; i < n; i++) if (fabs(A(i,i)) < tol) {
        printf("The element A(%i,%i)=%lf",i,i,A(i,i));
        return -2;
    }

    for (int k = 0; k < stepmax; k++) {
        maxerror = 0;
        for (int i = 0; i < n; i++) {
            sum = 0;
            for (int j = 0; j < n; j++) {
                if (i == j) continue;
                sum -= A(i,j) * x[j];
            }
            error = x[i];
            x[i] = 1.0/A(i,i) * ( b[i] + sum );
            error -= x[i];
            error = fabs(error);
            if (error > maxerror) maxerror = error;
            
        }

        if (maxerror < tol) return 0;

        if (k == stepmax-1) return -10;
        
    }

    #undef A

    return 0;
}