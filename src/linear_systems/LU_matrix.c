/*

This method tries to decompose an invertible matrix A into a product of two diagonal matrices LU,
Where L is a lower triangular matrix, and U is a upper triangular matrix, both non-singular.

It is known that every matrix can be decomposed to PA = T, where T echelon matrix and P represents
the transformations of the (arbitrary) matrix A. If we apply this result to a non-singular matrix A,
we have that P is a lower triangular matrix since it represents the operations from the highest row to the
lowest, and since those operations are exclusively applied to the lower rows, the matrix L has to be a 
lower triangular matrix. Since is a combination of elemental row operations, it is non-singular, and
the inverse is still a lower triangular matrix that represents the inverse operations, which are still
in the zone under the diagonal (changing the sign), so A = P^-1 * T, where P^-1 = L and T = U.

Of course, everything was assuming that the diagonal elements of A were different from zero. 
If any element were zero, we can exchange rows (another elemental operation) to make all of the
diagonal elements different from zero. The combination of row exchanges made (before substracting any row)
is called P, and so PA = LU for a permutation matrix P.

This is not a formal demostration, but it's a small approach of the meaning of the method.

Since PA = LU, then Pb = PAx = LUx, and the objective is to solve LUx = Pb, where L are the matrices
defined previously. Since both L and U are triangular matrices, with the change of variable Ux = y,
we can solve easily both systems, first Ly = Pb, and then Ux = y.

*/

#include <math.h>

int LU_matrix (int n, double *A, int *p, double tol) {

    #define A(i,j) A[i*n + j]

    int max;
    int temp;
    double m;

    // Computation of the permutation matrix (vector)

    for (int i = 0; i < n; i++) p[i] = i;

    for (int i = 0; i < n; i++){
        max = i;
        for (int j = i; j < n; j++){
            if(fabs(A(p[j],i)) > fabs(A(p[max],i))) max = j;
        }

        if (fabs(A(p[max],i)) < tol) return -1;

        temp = p[i];
        p[i] = p[max];
        p[max] = temp;

        for (int j = i+1; j < n; j++) {
            m = A(p[j],i)/A(p[i],i);
            for (int k = i+1; k < n; k++) A(p[j],k) -= m*A(p[i],k);
            A(p[j],i) = m;
        }

    }

    return 0;

}