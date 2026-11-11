#include <stdio.h>
#include "flow.h"
#include "shooting.h"

int shooting_function (int m, double t, double s[], double T[], void *prms) {

    shooting_prms *p = (shooting_prms *)prms;
    double x[p->n];
    double h = p->h;
    int (*vector_field) (int n, double t, double x[], double f[], void * prms) = p->vector_field;
    t = p->a;

    // s - initial values (known). Those are the values of the iteration
    // If we don't know the final time, s[n-1] is the final time 
    // which is also integrated.

    // x - the value of the flow with initial values s at time 1.

    // The shooting function is from Rm to Rm, where m is the number
    // of parameters we don't know.

    initial_shooting_conditions(x,s);

    flow(&t, x, &h, p->b, p->hmin, p->hmax, p->tol, p->stepmaxrk, p->n, vector_field, prms);


    // Write here the rest of the conditions

    T[0] = x[0] - 1;

    return 0;

}

int initial_shooting_conditions (double *x, double *s) {

    // Write here the initial conditions

    x[0] = 0;
    x[1] = s[0];

    return 0;

}