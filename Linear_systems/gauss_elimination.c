#include <stdio.h>
#include <math.h>
#include "matrix_print.h"

int gauss_elimination (int n, double *A, double *x, int pivot, double tol){

    #define A(i,j) A[i*(n+1)+j]

    double m = 0;
    double temp = 0;
    int maxrow = 0;
    int maxcolumn = 0;
    int permutation[n];

    for (int i=0;i<n;i++) permutation[i]=i;

    for (int i=0;i<n;i++){

        // Zero pivoting

        if (pivot==0 && i!=n-1) {

            if(fabs(A(i,i))<=tol) {
                for (int j=i+1;j<n;j++) {
                    if (fabs(A(j,i)) > tol) {
                        for (int k=i;k<n+1;k++) {
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

        if (pivot==1 && i!=n-1) {
            
            maxrow = i;
            for (int j=i+1;j<n;j++) {
                if (fabs(A(j,i)) > fabs(A(maxrow,i))) {
                    maxrow = j;
                }
            }

            for (int j=i;j<n+1;j++) {
                temp = A(i,j);
                A(i,j) = A(maxrow,j);
                A(maxrow,j) = temp;
            }

        } 

        // Complete pivoting

        if (pivot==2 && i!=n-1) {
            maxrow = i;
            maxcolumn = i;

            for (int j=i;j<n;j++){
                for (int k=i;k<n;k++){
                    if (fabs(A(j,permutation[k])) > fabs(A(maxrow,permutation[maxcolumn]))) {
                        maxrow = j;
                        maxcolumn = k;
                    }
                }
                
            }

            // Exchange rows

            for (int j=i;j<n;j++){
                temp = A(i,permutation[j]);
                A(i,permutation[j]) = A(maxrow,permutation[j]);
                A(maxrow,permutation[j]) = temp;
            }
            
            temp = A(i,n);
            A(i,n) = A(maxrow,n);
            A(maxrow,n) = temp;

            // Exchange columns

            temp = permutation[i];
            permutation[i] = permutation[maxcolumn];
            permutation[maxcolumn] = temp;

        }
        
        

        // Return error (-1) if it's still 0 (underdetermined/uncompatible linear system)

        if (fabs(A(i,permutation[i])) <= tol) return -1;

        for (int j=i+1;j<n;j++){

            m = A(j,permutation[i])/A(i,permutation[i]);

            for (int k=n-1;k>=i;k--) A(j,permutation[k]) -= m*A(i,permutation[k]);
            A(j,n) -= m*A(i,n);
        }
    }


    
    for (int i=n-1;i>=0;i--) {
        temp = A(i,n);
        for (int j=n-1;j>i;j--) {
            temp -= A(i,permutation[j])*x[permutation[j]];
        }
        x[permutation[i]]=temp/A(i,permutation[i]);
    }

    return 0;

}