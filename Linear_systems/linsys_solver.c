#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "matrix_print.h"
#include "gauss_elimination.h"
#include "LU_matrix.h"
#include "LU_solver.h"
#include "QR_matrix.h"
#include "QR_solver.h"

int main (int argc, char * argv[]) {

    int n, method, error;
    int *p;
    double tol;

    time_t start_time, finish_time;

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
    printf("4) LU decomposition\n");
    printf("5) QR decomposition\n");

    scanf("%i",&method);

    if (method <= 0 || method > 5) {
        printf("ERROR: Not a valid option\n");
        return 0;
    }


    FILE * system = fopen("system.txt","r");

    if (system == NULL) {
        printf("ERROR: System file couldn't be opened\n");
        return 0;
    }

    double *A;
    double *x;
    double *b;

    A = malloc((n*n)*sizeof(double));
    x = malloc(n*sizeof(double));
    b = malloc(n*sizeof(double));
    p = malloc(n*sizeof(int));
    
    for(int i = 0; i < n*n + n; i++){
        for (int j = 0; j < n; j++) fscanf(system, "%lf", &A[i*n + j]);
        fscanf(system, "%lf", &b[i]);
    }

    start_time = clock();

    if (method < 4) {
        error = gauss_elimination(n,A,x,b,(method-1)%3,tol);
    }

    if (method == 4) {
        error = LU_matrix(n,A,p,tol);
        if (error == -1) {
            printf("Underdetermined/incompatible system\n");
            return 0;
        }
        LU_solver(n,A,p,x,b,tol);
    }

    if (method == 5) {
        QR_matrix(n,n,A,p,tol);
        error = QR_solver(n,n,A,x,b,p,tol);
    }

    finish_time = clock();
    

    if (error == -1) {
        printf("Underdetermined/uncompatible system\n");
        return 0;

    }

    for (int i = 0; i < n; i++) printf("%.16g ",x[i]);
    printf("\n");

    printf("The algorithm has taken %.16g seconds\n",(double)(finish_time-start_time)/CLOCKS_PER_SEC);


    free(A);
    free(x);
    free(b);
    free(p);

    fclose(system);

    return 0;
}