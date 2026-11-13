/*

The shooting method is used mainly to solve numerically boundary value problems
(BVP) of ordinary differential equations. A boundary value problem is defined by:

x'(t) = F(t,x(t))
g(x(a),x(b)) = 0

where a,b are different real numbers (which usually represents a starting and
ending point), F : R^n --> R^n is the tangent field, and g : R^n x R^n --> R^n
is the function which defines the boundary condition.

The shooting method is based on the following logic: Consider an initial condition
x(a) = s where s is a real parameter, and this defines an unique solution x(t;s)
satisfying x(a;s) = s. Since the objetive is to satisfy g(x(a),x(b)) = 0,
it is equivalent to find s such that g(x(a),x(b)) = g(s,x(b;s)) = 0.
Then, we apply the Newton method to find the parameter s satisfying the boundary
condition.

In summary, we choose a seed s_0, we compute the IVP for x(a) = s_0, and the Newton 
method gives us another s_1 correcting the trajectory to minimize g(s,x(b;s)),
and do the same for x(a) = s_1.

It is extremely important to note that the existence nor uniqueness of such s is not guaranteed.

There are two ways of computing the jacobian of g: by variational equations or finite differences.
The election of either way is conditioned in the same way as the Newton method.
Variational equations ensures a more accuracy and a quadratic convergence, but derivatives
can be difficult depending on the problem. Finite differences are easier to compute, but we
lose a degree of convergence, being just linear.

*/

#include <stdio.h>
#include "newton_raphson.h"
#include "shooting.h"

int shooting_method (int m, double *hd, double *s, double *f, double *df, double tolsys, double tolnt,
    int (*shooting_function) (int n, double t, double x[], double f[], void *prms), 
    int (*shooting_jacobian) (int m, double t, double s[], double dT[], void * prms),
    void *prms) {

    // T(s) = g(s,x(b;s))
    
    return newton_raphson(m, hd, s, f, df, tolsys, tolnt, shooting_function, shooting_jacobian, prms);
    
}