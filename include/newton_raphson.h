int newton_raphson (int n, double *h, double x[], double f[], double df[], double tolsys, double tolnt,
    int (*vector_field) (int n, double t, double x[], double f[], void *prms),
    int (*jacobian) (int n, double t, double x[], double df[], void *prms),
    void * prms);