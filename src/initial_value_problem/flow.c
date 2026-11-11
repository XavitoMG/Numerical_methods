#include <stdio.h>
#include <math.h>
#include "rk78.h"

int flow (
    double *t, double x[], double *h, double T,
    double pasmin, double pasmax, double tol, int npasmx,
    int n, int (*camp)(int n, double t, double x[], double f[], void *prm),
    void *prm) 
{
    (*h) *= (T>=0?1:-1);
    double t0 = (*t);
    int err;

    for (int i = 0; fabs((*t) + (*h) - (t0)) <= fabs(T); i++) {

        err = rk78(t, x, h, pasmin, pasmax, tol, n, camp, prm);

        if (i == npasmx-1 || err < 0) {
            if (err < 0) printf("ERROR rk78\n");
            else printf("Step max reached\n");
            return -1;
        }
    }

    (*h) = T - (*t);
    err = rk78(t, x, h, pasmin, pasmax, tol, n, camp, prm);
    if (err < 0) {
        printf("ERROR rk78\n");
        return -1;
    }

    return 0;
      
}
