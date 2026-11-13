#include <stdio.h>
#include "flow.h"
#include "shooting.h"

int shooting_variationals (int dim, double t, double x[], double F[], void * prms);

int shooting_jacobian (int m, double t, double s[], double dT[], void * prms) {

    #define JT(i,j) dT[i*m+j]

    shooting_prms *p = (shooting_prms *)prms;
    int n = p->n;
    int k = p->k;
    double x[n+n*m];
    int assign[m];
    double h = p->h;

    initial_shooting_conditions(x,s,assign,prms);

    #define icx(i,j) x[n+i*(m-k)+j]
    #define icp(i,j) x[n+n*(m-k)+i*k+j]

    for (int i = 0; i < n; i++) for (int j = 0; j < m-k; j++) icx(i,j) = (i == assign[j] ? 1 : 0);
    for (int i = 0; i < n; i++) for (int j = 0; j < k; j++) icp(i,j) = 0;

    flow(&t,x,&h,p->b,p->hmin,p->hmax,p->tol,p->stepmaxrk,n+n*m,&shooting_variationals, prms);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m-k; j++) {
            JT(i,j) = x[n+i*(m-k)+j];
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            JT(i,(m-k) + j) = x[n + n*(m-k) + i*k + j];
        }
    }

    JT(1,0) += -1; 

    #undef JT
    #undef icx
    #undef icp

    return 0;
}

int shooting_variationals (int dim, double t, double x[], double F[], void * prms) {

    // n - number of dimensions of the original system
    // m - dimension of s
    // k - number of parameters that are not in the initial conditions
    // Thus, m-p is the number of free parameters in the initial conditions

    // n equations from the original system
    // n*(m-k) equations to find the derivatives in respect of the initial conditions
    // n*k equations to find the derivatives in respect of the other parameters
    // n+n*m equations in total

    shooting_prms *p = (shooting_prms *)prms;
    int k = p->k;
    int n = p->n;
    int m = p->m;

    double dfx[n*n];
    double dfp[n*k];

    #define f(i) F[i]
    #define D2f(i,j) dfx[i*n+j]
    #define D3f(i,j) dfp[i*k+j]
    #define FX(i,j) F[n+i*(m-k)+j]
    #define FP(i,j) F[n+n*(m-k)+i*k+j] 
    #define x(i) x[i]
    #define xs(i,j) x[n+i*(m-k)+j]
    #define xp(i,j) x[n+n*(m-k)+i*k+j]

    /*
    
    The variational equations are

    (x1)'
    (..)   = f(x), x0 from the iteration
    (xn)
    
    (x1_s1 ... x1_s(m-k))'         (x1_s1 ... x1_s(m-k))
    (...   ...    ...   )  = D2f(x)( ...  ...   ...    ),          x0 = (e_s1 ... e_s(m-k))
    (xn_s1 ... xn_s(m-k))          (xn_s1 ... xn_s(m-k))

    (x1_p1 ... x1_pk)'         (x1_p1 ... x1_pk)         (x1_p1 ... x1_pk)
    ( ...  ...  ... )  = D2f(x)( ...  ...  ... ) + D3f(x)( ...  ...  ... ),      x0 = 0_nxk
    (xn_p1 ... xn_pk)          (xn_p1 ... xn_pk)         (xn_p1 ... xn_pk)
 
    */

    // Declaration of the original function

    f(0) = p->T * x[1];
    f(1) = p->T * 2 * (1 - x[0] * x[0]) * x[1] - p->T * x[0];

    // Derivation in respect of x

    D2f(0,0) = 0;                               D2f(0,1) = p->T;
    D2f(1,0) = -p->T * (4 * x[0] * x[1] + 1);   D2f(1,1) = 2 * p->T * (1 - x[0] * x[0]);

    // Derivation in respect of parameters

    D3f(0,0) = x[1];
    D3f(1,0) = 2*(1 - x[0]*x[0]) * x[1] - x[0];

    // Declaration of the whole system

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m-k; j++) {
            FX(i,j) = 0;
            for (int l = 0; l < n; l++) FX(i,j) += D2f(i,l) * xs(l,j);
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            FP(i,j) = D3f(i,j);
            for (int l = 0; l < n; l++) FP(i,j) += D2f(i,l) * xp(l,j);
        }
    }


    #undef f
    #undef D2f
    #undef D3f
    #undef FX
    #undef FP
    #undef x
    #undef xs
    #undef xp

    return 0;
}