
int finite_difference (int n, double *h, double x[], double df[], 
    int (*vector_field) (int n, double t, double x[], double f[], void * prms), void * prms) {

    double f[n];
    double df_col[n];

    vector_field(n,0.0,x,f,prms);

    for (int i = 0; i < n; i++) {
        x[i] += *h;
        vector_field(n,0.0,x,df_col,prms);
        x[i] -= *h;
        for (int j = 0; j < n; j++) df[j*n+i] = (df_col[j] - f[j])/ *h;
    }

    return 0;
}