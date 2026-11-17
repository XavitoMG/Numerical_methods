#include <stdio.h>
#include <math.h>

#define MU 1.215058560962404e-2

int main (int argc, char *argv[]) {

    int n = 5;
    double x[n];

    double r0, rf, T, a, v0;

    r0 = 0.0170137704586768;
    rf = 0.02023161866935808;
    T = 0.8523191889448239;
    a = 0.9653659215446515;
    v0 = 10.67764174790542;


    x[0] = MU + r0 * cos(a);
    x[1] = r0 * sin(a);
    x[2] = -v0 * sin(a);
    x[3] = v0 * cos(a);
    x[4] = T;
    
    printf("%.16g %.16g %.16g %.16g %.16g\n",x[0],x[1],x[2],x[3],x[4]);

    

    return 0;

}