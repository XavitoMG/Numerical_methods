#include <math.h>

int LU_matrix(int n, double *A, int *p, double tol) {

    #define A(i,j) A[i*n + j]

    int max;
    int temp;
    double m;

    // Computation of the permutation matrix (vector)

    for (int i = 0; i < n; i++) p[i] = i;

    for (int i = 0; i < n; i++){
        max = i;
        for (int j = i; j < n; j++){
            if(fabs(A(p[j],i)) > fabs(A(p[max],i))) max = j;
        }

        if (fabs(A(p[max],i)) < tol) return -1;

        temp = p[i];
        p[i] = p[max];
        p[max] = temp;

        for (int j = i+1; j < n; j++) {
            m = A(p[j],i)/A(p[i],i);
            for (int k = i+1; k < n; k++) A(p[j],k) -= m*A(p[i],k);
            A(p[j],i) = m;
        }

    }

    return 0;

}