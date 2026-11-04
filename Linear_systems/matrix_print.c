#include <stdio.h>

int matrix_print(int n, int m, double *A){

    #define A(i,j) A[i*m+j]
    
    for (int i=0;i<n;i++){
        for (int j=0;j<m;j++) {
            printf("%.16g ",A(i,j));
        }
        printf("\n");
    }
}