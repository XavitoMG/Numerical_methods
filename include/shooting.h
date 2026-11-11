typedef struct 
{
    int n;
    double h;
    double hmin;
    double hmax;
    double tol;
    int stepmaxrk;
    double a;
    double b;
    int (*vector_field) (int n, double t, double *x, double *f, void *prms);
} shooting_prms;

int shooting_method (int m, double *hd, double *s, double *f, double *df, double tolsys, double tolnt,
    int (*shooting_function) (int n, double t, double x[], double f[], void *prms), 
    void *prms);

int shooting_function (int m, double t, double s[], double T[], void *prms);

int initial_shooting_conditions (double *x, double *s);