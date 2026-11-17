#include <math.h>
#include "shooting.h"

#define MU 1.215058560962404e-2

int vector_field (int n, double t, double x[], double f[], void *prms) {

    shooting_prms * p = (shooting_prms *)prms;

    double xmmu, xmmup1, y2, td1, td2, rh12, rh22, rh1, rh2,
    rh13, rh23, pot3, potx3, onemmu;
    onemmu=1-MU;
    xmmu=x[0]-MU; xmmup1=xmmu+1;
    y2=x[1]*x[1];
    rh12=xmmu*xmmu+y2; rh22=xmmup1*xmmup1+y2;
    rh1=sqrt(rh12); rh2=sqrt(rh22);
    rh13=rh12*rh1; rh23=rh22*rh2;
    td1=-onemmu/rh13; td2=-MU/rh23;
    pot3=td1+td2;
    potx3=td1*xmmu+td2*xmmup1;
    f[0]=p->T*x[2]; f[1]=p->T*x[3];
    f[2]=p->T*(2*x[3]+x[0]+potx3);
    f[3]=p->T*(-2*x[2]+x[1]*(1+pot3));

    return 0;
    
}