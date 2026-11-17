#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "shooting.h"
#include "vector_field.h"
#include "rk78.h"

#define MU 1.215058560962404e-2

int main (int argc, char * argv[]) {


    // WARNING: since the newton method is currently only programmed to solve
    // for square jacobians, if n!=m, the method won't work. 

    FILE *cts = fopen("data/constants_bvp.txt","r");
    FILE *sol = fopen("data/solution.txt","w");
    FILE *trj = fopen("data/trajectory.txt","w");

    shooting_prms prms;
    int n, m, stepmaxsh, conv;
    double hd, tolsys, tolsh, r0, rf, T, a, v0;
    double *f, *df;
    fscanf(cts,"%i %i %i",&n,&m,&prms.k);
    double s[m], x[n];

    f = malloc(m*sizeof(double));
    df = malloc(n*m*sizeof(double));

    fscanf(cts,"%lf",&prms.a);
    fscanf(cts,"%lf",&prms.b);
    for (int i = 0; i < m; i++) fscanf(cts,"%lf",s+i);
    fscanf(cts,"%lf",&prms.h);
    fscanf(cts,"%lf",&prms.hmin);
    fscanf(cts,"%lf",&prms.hmax);
    fscanf(cts,"%lf",&hd);
    fscanf(cts,"%lf",&tolsys);
    fscanf(cts,"%lf",&prms.tol);
    fscanf(cts,"%lf",&tolsh);
    fscanf(cts,"%i",&prms.stepmaxrk);
    fscanf(cts,"%i",&stepmaxsh);
    prms.n = n;
    prms.m = m;
    prms.vector_field = &vector_field_shooting;

    r0 = s[0];
    rf = s[1];
    T = s[2];
    a = s[3];
    v0 = s[4];

    prms.rf = rf;
    prms.r0 = r0;

    s[0] = MU + r0 * cos(a);
    s[1] = r0 * sin(a);
    s[2] = -v0 * sin(a);
    s[3] = v0 * cos(a);
    s[4] = T;

    for (int i = 0; i < m; i++) fprintf(sol, "%lf ", s[i]);
    fprintf(sol, "\n");

    for (int k = 0; k < stepmaxsh; k++) {

        printf("%.16g %.16g %.16g %.16g %.16g\n", s[0],s[1],s[2],s[3],s[4]);
        
        conv = shooting_method(m, &hd, s, f, df, tolsys, tolsh, &shooting_function, &shooting_jacobian, &prms);
        
        for (int i = 0; i < m; i++) fprintf(sol, "%.16g ", s[i]);
        fprintf(sol, "\n");

        if (conv == 0) {
            printf("The method has converged in %i iterations\n",k);
            break;
        }


        if (k == stepmaxsh-1) {
            printf("The method hasn't converged\n");
            return -1;
        }

    }

    initial_shooting_conditions(x,s,NULL,&prms);

    double t, h;
    t = prms.a;
    h = prms.h * (prms.b >= 0 ? 1 : -1);

    for (int i = 0; fabs(t - prms.a + h) < fabs(prms.b); i++){
        rk78(&t, x, &h, prms.hmin, prms.hmax, prms.tol, n, &vector_field, &prms);
        for (int j=0; j<n; j++) fprintf(trj, "%.16g ", x[j]);
        fprintf(trj, "\n");
    }

    h = prms.b - t;
    rk78(&t, x, &h, prms.hmin, prms.hmax, prms.tol, n, &vector_field, &prms);
    for (int j=0; j<n; j++) fprintf(trj, "%.16g ", x[j]);


    return 0;

    
} 