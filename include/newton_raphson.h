int newton_raphson (int n, double *h, double x[], double f[], double df[], double tolsys, 
    int (*VF) (int n, double t, double x[], double f[], void * prms), void * prms);