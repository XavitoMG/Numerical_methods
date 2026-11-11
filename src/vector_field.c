#include <math.h>

int vector_field (int n, double t, double x[], double f[], void * prms) {

    f[0] = x[1];
    f[1] = - x[0] * x[0] * x[0];

    return 0;

}