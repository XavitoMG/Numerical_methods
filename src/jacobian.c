#include <math.h>

int jacobian (int n, double t, double x[], double df[], void * prms) {

    #define df(i,j) df[i*n + j]

    df(0,0) = 4*pow(x[0],4) - 6 * x[0] * x[1];        df(0,1) = -3 * pow(x[1],2) + 4 * pow(x[1],3);
    df(1,0) = 1/x[0] + 2 * pow(x[0],2) * pow(x[1],2); df(1,1) = 2 * pow(x[0],3) * x[1];

    #undef df

    return 0;
}