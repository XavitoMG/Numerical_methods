#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rk78.h"
#include "camp.h"

int main(int argc, char *argv[]) {

    FILE *cts = fopen("constants.txt","r");
    FILE *trj = fopen("trajectoria2.txt","w");

    int n, pn;
    double t0, t, h, hmin, hmax, tol, T;

    fscanf(cts, "%i %lf %lf", &n, &t0, &T);
    
    double x[n];

    for (int i=0; i<n; i++) fscanf(cts, "%lf", &x[i]);
    fscanf(cts, "%lf %lf %lf %lf", &h, &hmin, &hmax, &tol);
    fscanf(cts, "%i", &pn);

    double * p = (double *) malloc(pn * sizeof(double));

    for (int i=0; i<pn; i++) fscanf(cts, "%lf", &p[i]);

    t = t0;
    h *= (T>=0?1:-1);

    for (int i=0; fabs(t-t0+h)<fabs(T); i++){
        rk78(&t, x, &h, hmin, hmax, tol, n, &CV, p);
        for (int j=0; j<n; j++) fprintf(trj, "%.16g ", x[j]);
        fprintf(trj, "\n");
    }

    h=T-t;
    rk78(&t, x, &h, hmin, hmax, tol, n, &CV, p);
    for (int j=0; j<n; j++) fprintf(trj, "%.16g ", x[j]);
    
    

}