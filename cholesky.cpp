#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include "../linearalgebra.h"

/* -------------------------- cholesky -------------------------- */
/*  Given a hermitian positive definite matrix and its dimension,
    this function computes an upper triangular matrix R such that 
    A = R^* R.

    Input variables:
        a : pointer to array whose elements are the columns
             of hermitian positive definite matrix A.
        m : number of rows and number of columns in matrix A.

    Features: This algorithm has time complexity ~(1/3)m^3
    flops and requires O(1) additional memory.                   */

void cholesky (double ** a, int m) {
    int i, j;
    double t;
    
    for(i = 0; i < m; i++) {
        for(j = i + 1; j < m; j++) {
            /* a[j:m][j] -= a[j:m][i] a[i][j] / a[i][i] */
            rowsubrow(a, a[i][j] / a[i][i], m - j, j, j, i);   
        }
        /* a[i:m][i] = a[i:m][i]/sqrt(a[i][i]) */
        t = sqrt(a[i][i]);
        subscalar_div(a[i], t, i, m);
    }
}

int main () {
    int i, j, m;
    double x;

#if 0
    /* let user set the dimension of matrix A */
    printf("Enter the dimension m (where A is a m by m matrix): ");
    scanf("%i", &m);
#endif

    m = 3; 

    /* allocate memory for A */
    double ** a = new double*[m];
    for(i = 0; i < m; i++) {
        a[i] = new double[m];
    }

    /* initialize the values in matrix A */
    a[0][0] = 2;    a[1][0] = -1;   a[2][0] = 0;
    a[0][1] = -1;   a[1][1] = 2;    a[2][1] = -1;
    a[0][2] = 0;    a[1][2] = -1;   a[2][2] = 2;
    
#if 0
    for(i = 0; i < m; i++) {
        for(j = 0; j < m; j++) {
            if(j < i) {
                a[i][j] = i - j + 1; // this choice of values was arbitrary
            }
            else {
                a[i][j] = j - i + 1; // chosen so symmetric matrix
            }
        }
    }
#endif

    /* print the matrix A before calling cholesky */
    printf("A = \n");
    for(i = 0; i < m; i++) {
        for(j = 0; j < m; j++) {

            printf("%9.6g ", a[j][i]);
        }
        printf("\n");
    }
    printf("\n");

    /* execute cholesky */
    cholesky(a, m);

    /* print the matrix R (stored in A) after calling cholesky */
    printf("R* = \n");
    for(i = 0; i < m; i++) {
        for(j = 0; j < m; j++) {
            if(j < i) {
                printf("%9.6g ", 0.0);
            }
            else{
                printf("%9.6g ", a[i][j]);
            }
        }
        printf("\n");
    }
    printf("\n");

    /* free memory */
    for(i = 0; i < m; i++) {
        delete[] a[i];
    }
    delete[] a;
    return 0;
}
