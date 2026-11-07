#include "LU_decomposition.h"
#include "derivatives.h"

int newton_raphson (int n, double *h, double x[], double f[], double df[], double tolsys,
    int (*VF) (int n, double x[], double f[], void * prms), void * prms) {

    int p[n];
    double u[n];

    VF(n,x,f,prms);

    for (int i = 0; i < n; i++) f[i] *= -1.0;

    if ( *h > 0 ) finite_difference(n,h,x,df,VF,prms);

    if ( *h == 0) eval_derivative(n,x,df,prms);

    LU_matrix(n,df,p,tolsys);

    LU_solver(n,df,p,u,f,tolsys);

    for (int i = 0; i < n; i++) x[i] += u[i];
    
    return 0;

}