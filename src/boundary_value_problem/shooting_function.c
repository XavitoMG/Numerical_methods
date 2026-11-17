#include <stdio.h>
#include <math.h>
#include "flow.h"
#include "shooting.h"

#define MU 1.215058560962404e-2

int vector_field_shooting (int n, double t, double x[], double f[], void *prms) {

    shooting_prms * p = (shooting_prms *)prms;

    double xmmu, xmmup1, y2, td1, td2, rh12, rh22, rh1, rh2,
    rh13, rh23, pot3, potx3, onemmu;
    onemmu=1-MU;
    xmmu=x[0]-MU; xmmup1=xmmu+1;
    y2=x[1]*x[1];
    rh12=xmmu*xmmu+y2; rh22=xmmup1*xmmup1+y2;
    rh1=sqrt(rh12); rh2=sqrt(rh22);
    rh13=rh12*rh1; rh23=rh22*rh2;
    td1=-onemmu/rh13; td2=-MU/rh23;
    pot3=td1+td2;
    potx3=td1*xmmu+td2*xmmup1;
    f[0]=p->T*x[2]; f[1]=p->T*x[3];
    f[2]=p->T*(2*x[3]+x[0]+potx3);
    f[3]=p->T*(-2*x[2]+x[1]*(1+pot3));

    return 0;
    
}

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

    initial_shooting_conditions(x,s,NULL,prms);

    flow(&t, x, &h, p->b, p->hmin, p->hmax, p->tol, p->stepmaxrk, p->n, vector_field, prms);


    // Write here the rest of the conditions

    
    T[0] = (s[0] - MU) * (s[0] - MU) + s[1] * s[1] - p->r0 * p->r0;
    T[1] = (s[0] - MU) * s[2] + s[1] * s[3];
    T[2] = x[1];
    T[3] = x[2];
    T[4] = x[0] - (MU - 1 - p->rf);

    return 0;

}

int initial_shooting_conditions (double *x, double *s, int *assign, void * prms) {

    // Write here the initial conditions

    shooting_prms *p = (shooting_prms *)prms;

    x[0] = s[0];
    x[1] = s[1];
    x[2] = s[2];
    x[3] = s[3];
    p->T = s[4];

    if (assign != NULL) {
        assign[0] = 0;
        assign[1] = 1;
        assign[2] = 2;
        assign[3] = 3;
    }
    

    return 0;

}