#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "vector_field.h"
#include "newton_raphson.h"

int main (int argc, char * argv[]) {

    if (argc != 6) {
        printf("/main n h tolsys tolnt stepmax\n");
        return 0;
    }

    int n = atoi(argv[1]);
    double h = atof(argv[2]); 
    double tolsys = atof(argv[3]);
    double tolnt = atof(argv[4]);
    int stepmax = atoi(argv[5]);

    double xprev[n], x[n], f[n], *df, errmaxx, errmaxf, errx, errf;

    FILE * seed = fopen("seed.txt","r");

    FILE * solution = fopen("solution.txt","w");

    df = malloc(n*n*sizeof(double));

    for (int i = 0; i < n; i++) fscanf(seed,"%lf",&x[i]);


    for (int i = 0; i < n; i++) fprintf(solution,"%.16g ",x[i]);
    fprintf(solution,"\n");

    for (int k = 0; k < stepmax; k++) {

        memcpy(xprev,x,n*sizeof(double));
        newton_raphson(n,&h,x,f,df,tolsys,&vector_field,NULL);

        for (int i = 0; i < n; i++) fprintf(solution,"%.16g ",x[i]);
        fprintf(solution,"\n");
        
        errmaxx = 0;
        errmaxf = 0;
        for (int i = 0; i < n; i++) {
            errx = fabs(x[i] - xprev[i]);
            errf = fabs(f[i]);
            if (errx > errmaxx) errmaxx = errx;
            if (errf > errmaxf) errmaxf = errf;
        }

        if (errmaxx < tolnt || errmaxf < tolnt) {
            printf("The method has converged in %i iterations\n",k+1);
            return 0;
        }
        
        if (k == stepmax - 1) {
            printf("The method did not converge in %i iterations\n",stepmax);
            return 0;
        }
    }

}