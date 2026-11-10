#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "matrix_print.h"
#include "linsyslib.h"

int main (int argc, char * argv[]) {

    int n, method, error, stepmax;
    int *p;
    double tol;

    time_t start_time, finish_time;

    if (argc != 4) {

        printf("./linsys_solver n tol stepmax\n");

        return 0;
    }

    tol = atof(argv[2]);
    n = atoi(argv[1]);
    stepmax = atoi(argv[3]);
    
    printf("Choose the numerical method to solve the linear system:\n");
    printf("1) Gauss elimination without pivotage\n");
    printf("2) Gauss elimination with maximum pivotage\n");
    printf("3) Gauss elimination with full pivotage\n");
    printf("4) LU decomposition\n");
    printf("5) QR decomposition\n");
    printf("6) Gauss-Seidel method (iterative)\n");

    scanf("%i",&method);

    if (method <= 0 || method > 6) {
        printf("ERROR: Not a valid option\n");
        return 0;
    }


    FILE * system = fopen("data/system.txt","r");

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
    
    for(int i = 0; i < n; i++){
        for (int j = 0; j < n; j++) fscanf(system, "%lf", &A[i*n + j]);
        fscanf(system, "%lf", &b[i]);
    }

    if (method > 4) {
        for (int i = 0; i < n; i++) fscanf(system, "%lf", x+i);
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

    if (method == 6) {
        error = gauss_seidel(n,A,x,b,tol,stepmax);
    }

    finish_time = clock();
    

    if (error == -1) {
        printf("Underdetermined/uncompatible system\n");
        return 0;
    }
    if (error == -2) {
        printf("An element from the diagonal is 0 (jacobi can't proceed)\n");
        return 0;
    }
    if (error == -10) {
        printf("The method did not converge\n");
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