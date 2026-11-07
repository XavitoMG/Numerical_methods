int gauss_elimination (int n, double *A, double *x, double *b, int pivot, double tol);

int LU_matrix (int n, double *A, int *p, double tol);
int LU_solver (int n, double *A, int *p, double *x, double *b, double tol);

int QR_matrix (int n, int m, double *A, int *p, double tol);
int QR_solver (int n, int m, double *A, double *x, double *b, int *p, double tol);

int gauss_seidel (int n, double *A, double *x, double *b, double tol, int stepmax);