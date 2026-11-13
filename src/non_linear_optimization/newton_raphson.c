/*

Consider a non-linear system of equations F(x)=0.
The Newton-Raphson method is based on the Taylor expansion of F, which is

F(x) = F(x_0) + JF(x_0)*(x - x_0) + O(x^2)

Evaluating F at a zero, x=a, we get

0 = F(a) = F(x) + JF(x)*(a - x)

If we consider x_n+1 an approximation of a, we could construct an iterative method to find a,

JF(x_n)*(x_n+1 - x_n) = -F(x_n)

With the change of variable u_n = x_n+1 - x_n, the problem reduces to solve the next linear system:

JF(x_n)u = -F(x_n)

and then reverse the change of variable u = u_n = x_n+1 - x_n.

The jacobian can be computed either by an approximation by finite difference method (linear convergence, easier to compute) or
by computing the matrix analytically (quadratic convergence, can be difficult or not worth it sometimes).

*/

#include <math.h>
#include <string.h>
#include "linsyslib.h"
#include "derivatives.h"
#include "matrix_print.h"

int newton_raphson (int n, double *h, double x[], double f[], double df[], double tolsys, double tolnt,
    int (*vector_field) (int n, double t, double x[], double f[], void *prms),
    int (*jacobian) (int n, double t, double x[], double df[], void *prms),
    void * prms)
    {

    int p[n];
    double u[n];

    vector_field(n,0.0,x,f,prms);

    double xprev[n], errmaxx, errmaxf, errx, errf;

    memcpy(xprev,x,n*sizeof(double));

    if ( *h > 0 ) finite_difference(n,h,x,f,df,vector_field,prms);

    if ( *h == 0) jacobian(n,0.0,x,df,prms);

    //matrix_print(n,n,df);

    for (int i = 0; i < n; i++) f[i] *= -1.0;

    LU_matrix(n,df,p,tolsys);

    LU_solver(n,df,p,u,f,tolsys);


    for (int i = 0; i < n; i++) x[i] += u[i];

    errmaxx = 0;
    errmaxf = 0;
    for (int i = 0; i < n; i++) {
        errx = fabs(x[i] - xprev[i]);
        errf = fabs(f[i]);
        if (errx > errmaxx) errmaxx = errx;
        if (errf > errmaxf) errmaxf = errf;
    }

    if (errmaxx < tolnt || errmaxf < tolnt) return 0;
    
    return -1;

}