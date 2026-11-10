int finite_difference (int n, double *h, double x[], double *df, 
    int (*vector_field) (int n, double t, double x[], double f[], void * prms), void * prms);
int eval_derivative (int n, double x[], double *df, void * prms);