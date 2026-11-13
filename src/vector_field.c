#include <math.h>
#include "shooting.h"

int vector_field (int n, double t, double x[], double f[], void * prms) {

    shooting_prms * p = (shooting_prms *)prms;

    f[0] = p->T * x[1];
    f[1] = 2*p->T*x[1]*(1-x[0]*x[0]) - p->T*x[0];

    return 0;

}