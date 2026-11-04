#include <stdio.h>
#include <stdlib.h>
#include "gauss_elimination.h"

int main(int argc, char * argv[]) {

    int n, method, error;
    double tol;

    if (argc != 3) {

        printf("./linsys_solver n tol\n");

        return 0;
    }

    tol = atof(argv[2]);
    n = atoi(argv[1]);
    
    printf("Choose the numerical method to solve the linear system:\n");
    printf("1) Gauss elimination without pivotage\n");
    printf("2) Gauss elimination with maximum pivotage\n");
    printf("3) Gauss elimination with full pivotage\n");

    scanf("%i",&method);

    if (method<=0 || method>3) {
        printf("ERROR: Not a valid option\n");
        return 0;
    }


    FILE * system = fopen("system.txt","r");

    double *A;
    double *x;

    A = malloc((n*n+n)*sizeof(double));
    x = malloc(n*sizeof(double));
    
    for(int i=0;i<n*n+n;i++){
        fscanf(system,"%lf",A+i);
    }

    if (method<4) {
        error = gauss_elimination(n,A,x,(method-1)%3,tol);
    }
    
    

    if (error == -1) {
        printf("Underdetermined/uncompatible system\n");
        return 0;

    }

    for (int i=0;i<n;i++) printf("%.16g ",x[i]);
    printf("\n");

}