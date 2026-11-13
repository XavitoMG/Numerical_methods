typedef struct 
{
    int n;
    int m;
    int k;
    double h;
    double hmin;
    double hmax;
    double tol;
    int stepmaxrk;
    double a;
    double b;
    int (*vector_field) (int n, double t, double *x, double *f, void *prms);
    double T;
} shooting_prms;

int shooting_method (int m, double *hd, double *s, double *f, double *df, double tolsys, double tolnt,
    int (*shooting_function) (int n, double t, double x[], double f[], void *prms), 
    int (*shooting_jacobian) (int m, double t, double s[], double dT[], void *prms),
    void *prms);

int vector_field_shooting (int n, double t, double x[], double f[], void *prms);

int shooting_function (int m, double t, double s[], double T[], void *prms);

int initial_shooting_conditions (double *x, double *s, int *assign, void *prms);

int shooting_jacobian (int m, double t, double s[], double dT[], void *prms);