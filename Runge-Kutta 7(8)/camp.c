struct CV_prms {
    double a,b,c;
};

int CV (int n, double t, double x[], double f[], void * prm) {

    struct CV_prms *cvp = (struct CV_prms *) prm;

    f[0] = cvp->a*(x[1] - x[0]);
    f[1] = x[0]*(cvp->b - x[2]) - x[1];
    f[2] = x[0]*x[1] - cvp->c*x[2];

    return 0;
}