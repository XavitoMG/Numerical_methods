#include <math.h>

int eval_derivative (int n, double x[], double df[], void * prms) {

    #define df(i,j) df[i*n + j]

    df(0,0) = exp(x[0]-1) + x[1];  df(0,1) = x[0]      ;  df(0,2) = 0;
    df(1,0) = 1                 ;  df(1,1) = cos(x[1]) ;  df(1,2) = 2*x[2];
    df(2,0) = 3 * x[0] * x[0]   ;  df(2,1) = -exp(x[1]);  df(2,2) = 1;

    #undef df

    return 0;
}