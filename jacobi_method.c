/*
-In numerical linear algebra, the Jacobi method determining the solutions of a diagonally dominant system of linear equations.
-This method is based on the transformation of the linear system of a Ax = b into the system x = Cx + d, in which the matrix C has zeros on the diagonal.
-The vector is updated using previous estimate for all components of x to evaluate the right hand side of the equation.

*/

#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <string.h>  

#define MAX_ITER 1000      // Maksimum iterasyon sayisi
#define TOLERANCE 1e-6     // Hata toleransi

int iterations = 0;

void jacobi_method(int n, double **a, double *b, double *x) {
    double *x_new = (double*)malloc(n * sizeof(double));
    double error;

    do {
        error = 0.0;
        for (int i = 0; i < n; i++) {
            double sum = b[i];
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    sum -= a[i][j] * x[j];
                }
            }
            x_new[i] = sum / a[i][i];
            error += fabs(x_new[i] - x[i]);
        }
        for (int i = 0; i < n; i++) {
            x[i] = x_new[i];
        }
        for (int i = 0; i < n; i++) {
        printf("x[%d] = %.5f\t", i, x[i]);
       }
       printf("\n");
        iterations++;
    } while (error > TOLERANCE && iterations < MAX_ITER);

    free(x_new);
}

int main() {
    int n;

    printf("Please enter matrix size(n): ");
    scanf("%d", &n);


    double **a = (double**)malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) {
        a[i] = (double*)malloc(n * sizeof(double));
    }
    double *b = (double*)malloc(n * sizeof(double));
    double *x = (double*)malloc(n * sizeof(double));

    // x dizisinin degerlerini sifirladik
    memset(x, 0, n * sizeof(double));

    // Katsayılar matrisini ve b vektörünü aldik
    printf("Please enter coefficient matrix (A)\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("A[%d][%d]: ", i, j);
            scanf("%lf", &a[i][j]);
        }
    }

    printf("Enter the right hand side of the equations (b):\n");
    for (int i = 0; i < n; i++) {
        printf("b[%d]: ", i);
        scanf("%lf", &b[i]);
    }

    // Jacobi metodunu cagirdik
    jacobi_method(n, a, b, x);

    // Çözümleri ekrana yazdirdik
    printf("\nSolutions (x):\n");
    for (int i = 0; i < n; i++) {
        printf("x[%d] = %.5f\n", i, x[i]);
    }

    printf("At iteration %d, tolerance is almost zero!", iterations);

    // Belleği serbest biraktik
    for (int i = 0; i < n; i++) {
        free(a[i]);
    }
    free(a);
    free(b);
    free(x);

    return 0;
}

