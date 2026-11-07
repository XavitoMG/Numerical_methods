#include <math.h>

int vector_field (int n, double x[], double f[], void * prms) {

    f[0] = exp(x[0] - 1) + x[0] * x[1] - 1;
    f[1] = sin(x[1]) + x[2] * x[2] + x[0] - 1;
    f[2] = x[0] * x[0] * x[0] + x[2] - exp(x[1]);

    return 0;

}