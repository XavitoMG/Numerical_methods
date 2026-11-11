#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "shooting.h"
#include "vector_field.h"
#include "rk78.h"

int main (int argc, char * argv[]) {


    FILE *cts = fopen("data/constants_bvp.txt","r");
    FILE *sol = fopen("data/solution.txt","w");
    FILE *trj = fopen("data/trajectory.txt","w");

    int n, m, pn, stepmaxsh, conv;
    double hd, tolsys, tolsh;
    double *f, *df;
    fscanf(cts,"%i %i %i",&pn,&n,&m);
    double s[m], x[n];

    shooting_prms prms;

    f = malloc(m*sizeof(double));
    df = malloc(m*m*sizeof(double));

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
    prms.vector_field = &vector_field;

    for (int i = 0; i < m; i++) fprintf(sol, "%lf ", s[i]);
    fprintf(sol, "\n");

    for (int k = 0; k < stepmaxsh; k++) {
        
        conv = shooting_method(m, &hd, s, f, df, tolsys, tolsh, &shooting_function, &prms);
        
        for (int i = 0; i < m; i++) fprintf(sol, "%lf ", s[i]);
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

    initial_shooting_conditions(x,s);

    double t, h;
    t = prms.a;
    h = prms.h * (prms.b >= 0 ? 1 : -1);

    for (int i = 0; fabs(t - prms.a + h) < fabs(prms.b); i++){
        rk78(&t, x, &h, prms.hmin, prms.hmax, prms.tol, n, &vector_field, NULL);
        for (int j=0; j<n; j++) fprintf(trj, "%.16g ", x[j]);
        fprintf(trj, "\n");
    }

    h = prms.b - t;
    rk78(&t, x, &h, prms.hmin, prms.hmax, prms.tol, n, &vector_field, NULL);
    for (int j=0; j<n; j++) fprintf(trj, "%.16g ", x[j]);


    return 0;

    
} 