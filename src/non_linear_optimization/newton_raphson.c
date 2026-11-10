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

#include "linsyslib.h"
#include "derivatives.h"

int newton_raphson (int n, double *h, double x[], double f[], double df[], double tolsys,
    int (*VF) (int n, double t, double x[], double f[], void * prms), void * prms) {

    int p[n];
    double u[n];

    VF(n,0.0,x,f,prms);

    for (int i = 0; i < n; i++) f[i] *= -1.0;

    if ( *h > 0 ) finite_difference(n,h,x,df,VF,prms);

    if ( *h == 0) eval_derivative(n,x,df,prms);

    LU_matrix(n,df,p,tolsys);

    LU_solver(n,df,p,u,f,tolsys);

    for (int i = 0; i < n; i++) x[i] += u[i];
    
    return 0;

}