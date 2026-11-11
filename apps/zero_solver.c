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

    double x[n], f[n], *df;
    int conv;

    FILE * seed = fopen("data/seed.txt","r");

    FILE * solution = fopen("data/solution.txt","w");

    df = malloc(n*n*sizeof(double));

    for (int i = 0; i < n; i++) fscanf(seed,"%lf",&x[i]);


    for (int i = 0; i < n; i++) fprintf(solution,"%.16g ",x[i]);
    fprintf(solution,"\n");

    for (int k = 0; k < stepmax; k++) {

        conv = newton_raphson(n,&h,x,f,df,tolsys,tolnt,&vector_field,NULL);

        for (int i = 0; i < n; i++) fprintf(solution,"%.16g ",x[i]);
        fprintf(solution,"\n");

        if (conv == 0) {
            printf("The method has converged in %i terations\n",k);
            return 0;
        }

        else if (k == stepmax-1) {
            printf("The method has not converged\n");
            return -1;
        }
    }

}