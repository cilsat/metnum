#include "num.h"
#include <stdio.h>

#define PREC "%.3f "

void reset(double *arr, long n) {
    for (long i = 0; i < n; i++) {
        printf(PREC, arr[i]);
        arr[i] = 0.f;
    }
}

int main() {
    long i, n;
    matrix *a;
    double *b, *c;

    /***** SOAL 1 *****/
    n = 4;
    a = m_init(n, n);
    b = (double *) malloc(n * sizeof(double));
    c = (double *) malloc(n * sizeof(double));

    for (i = 0; i < n; i++) {
        c[i] = 0.f;
    }

    a->data[0][0] = 0.31; a->data[0][1] = 0.14; a->data[0][2] = 0.30; a->data[0][3] = 0.27; a->data[1][0] = 0.26; a->data[1][1] = 0.32; a->data[1][2] = 0.18; a->data[1][3] = 0.24; a->data[2][0] = 0.61; a->data[2][1] = 0.22; a->data[2][2] = 0.20; a->data[2][3] = 0.31; a->data[3][0] = 0.4; a->data[3][1] = 0.34; a->data[3][2] = 0.36; a->data[3][3] = 0.17; 
    b[0] = 1.02;
    b[1] = 1.00;
    b[2] = 1.34;
    b[3] = 1.27;

    printf("\nsoal 1\n");
    m_print(a);
    printf("gauss naive:\n ");
    gaussnaive(a, b, c);
    reset(c, n);
    printf("\n");
    printf("gauss jordan:\n ");
    gaussjordan(a, b, c);
    reset(c, n);
    printf("\n");
    printf("crout decomposition:\n ");
    crout(a, b, c);
    reset(c, n);
    printf("\n");
    printf("doolittle decomposition:\n ");
    doolittle(a, b, c);
    reset(c, n);
    printf("\n");

    m_del(a);
    free(b);
    free(c);

    /***** SOAL 2 *****/
    n = 9;
    a = m_init(n, n);
    b = (double *) malloc(n * sizeof(double));
    c = (double *) malloc(n * sizeof(double));

    for (i = 0; i < n; i++) {
        c[i] = 0.f;
    }

    a->data[0][0] = 0.7071; 
    a->data[0][3] = -1;
    a->data[0][3] = -0.8660;
    a->data[1][0] = 0.7071;
    a->data[1][2] = 1;
    a->data[1][4] = 0.5;
    a->data[2][1] = 1;
    a->data[2][5] = -1;
    a->data[3][2] = -1;
    a->data[4][6] = 1;
    a->data[4][8] = 0.7071;
    a->data[5][3] = 1;
    a->data[5][8] = -0.7071;
    a->data[6][4] = 0.8660;
    a->data[6][5] = 1;
    a->data[6][7] = -1;
    a->data[7][4] = -0.5;
    a->data[4][6] = -1;
    a->data[8][8] = 0.7071;
    b[1] = -1000;
    b[4] = 500;
    b[7] = -500;

    printf("\nsoal 2\n");
    m_print(a);
    printf("gauss naive:\n ");
    gaussnaive(a, b, c);
    reset(c, n);
    printf("\n");
    printf("gauss jordan:\n ");
    gaussjordan(a, b, c);
    reset(c, n);
    printf("\n");
    printf("doolittle decomposition:\n ");
    doolittle(a, b, c);
    reset(c, n);
    printf("\n");
    printf("crout decomposition:\n ");
    crout(a, b, c);
    reset(c, n);
    printf("\n");

    m_del(a);
    free(b);
    free(c);

    /**SOAL 3**/
    n = 6;
    a = m_init(n, n);
    b = (double *) malloc(n * sizeof(double));
    c = (double *) malloc(n * sizeof(double));

    for (i = 0; i < n; i++) {
        b[i] = 1.f;
        c[i] = 0.f;
    }

    printf("\nsoal 3\n");
    m_hilbert(a);
    m_print(a);
    printf("gauss naive:\n ");
    gaussnaive(a, b, c);
    reset(c, n);
    printf("\n");
    printf("gauss jordan:\n ");
    gaussjordan(a, b, c);
    reset(c, n);
    printf("\n");
    printf("doolittle decomposition:\n ");
    doolittle(a, b, c);
    reset(c, n);
    printf("\n");
    printf("crout decomposition:\n ");
    crout(a, b, c);
    reset(c, n);
    printf("\n");

    m_del(a);
    free(b);
    free(c);

    /** SOAL 4 (TERAPAN) **/
    n = 6;
    a = m_init(n, n);
    b = (double *) malloc(n * sizeof(double));
    c = (double *) malloc(n * sizeof(double));

    for (i = 0; i < n; i++) {
        c[i] = 0.f;
    }

    a->data[0][0] = 0.866;
    a->data[0][2] = -0.5;
    a->data[1][0] = 0.5;
    a->data[1][2] = 0.866;
    a->data[2][0] = -0.866;
    a->data[2][1] = -1;
    a->data[2][3] = -1;
    a->data[3][0] = -0.5;
    a->data[3][4] = -1;
    a->data[4][1] = 1;
    a->data[4][2] = 0.5;
    a->data[5][2] = -0.866;
    a->data[5][5] = -1;
    b[1] = -1000;

    printf("\nsoal terapan\n");
    m_print(a);
    printf("doolittle decomposition:\n ");
    doolittle(a, b, c);
    reset(c, n);
    printf("\n");
    printf("crout decomposition:\n ");
    crout(a, b, c);
    reset(c, n);
    printf("\n");

    m_del(a);
    free(b);
    free(c);

    return 0;
}
