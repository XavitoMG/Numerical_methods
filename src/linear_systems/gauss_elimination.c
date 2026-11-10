#include <stdio.h>
#include <math.h>
#include "matrix_print.h"

int gauss_elimination (int n, double *A, double *x, double *b, int pivot, double tol){

    #define A(i,j) A[i*n + j]

    double m = 0;
    double temp = 0;
    int maxrow = 0;
    int maxcolumn = 0;
    int p[n];

    for (int i = 0; i < n; i++) p[i]=i;

    for (int i = 0; i < n; i++){

        // Zero pivoting

        if (pivot == 0 && i != n-1) {

            if (fabs(A(i,i)) <= tol) {
                for (int j = i+1; j < n; j++) {
                    if (fabs(A(j,i)) > tol) {
                        for (int k = i; k < n + 1; k++) {
                            temp = A(i,k);
                            A(i,k) = A(j,k);
                            A(j,k) = temp;
                        }
                        break;
                    }
                }
            }
        }

        // Maximum pivoting

        if (pivot == 1 && i != n-1) {
            
            maxrow = i;
            for (int j = i+1; j < n; j++) {
                if (fabs(A(j,i)) > fabs(A(maxrow,i))) {
                    maxrow = j;
                }
            }

            for (int j = i; j < n; j++) {
                temp = A(i,j);
                A(i,j) = A(maxrow,j);
                A(maxrow,j) = temp;
            }

            temp = b[i];
            b[i] = b[maxrow];
            b[maxrow] = temp;


        } 

        // Complete pivoting

        if (pivot==2 && i!=n-1) {
            maxrow = i;
            maxcolumn = i;

            for (int j = i; j < n; j++){
                for (int k = i; k < n; k++){
                    if (fabs(A(j,p[k])) > fabs(A(maxrow,p[maxcolumn]))) {
                        maxrow = j;
                        maxcolumn = k;
                    }
                }
                
            }

            // Exchange rows

            for (int j = i; j < n; j++){
                temp = A(i,p[j]);
                A(i,p[j]) = A(maxrow,p[j]);
                A(maxrow,p[j]) = temp;
            }
            
            temp = b[i];
            b[i] = b[maxrow];
            b[maxrow] = temp;

            // Exchange columns

            temp = p[i];
            p[i] = p[maxcolumn];
            p[maxcolumn] = temp;

        }
        
        

        // Return error (-1) if it's still 0 (underdetermined/uncompatible linear system)

        if (fabs(A(i,p[i])) <= tol) return -1;

        for (int j = i+1; j < n; j++){

            m = A(j,p[i])/A(i,p[i]);

            for (int k = n-1; k >= i; k--) A(j,p[k]) -= m*A(i,p[k]);

            b[j] -= m*b[i];

        }
    }

    


    
    for (int i=n-1;i>=0;i--) {
        temp = b[i];
        for (int j=n-1;j>i;j--) temp -= A(i,p[j])*x[p[j]];
        x[p[i]]=temp/A(i,p[i]);
    }

    return 0;

}