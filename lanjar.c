#include "num.h"
#include <stdio.h>

int main() {
    long i, size;
    matrix *a;
    double *b, *c;

    /***** SOAL 1 *****/
    size = 4;
    a = m_init(size, size);
    b = (double *) malloc(size * sizeof(double));
    c = (double *) malloc(size * sizeof(double));

    for (i = 0; i < size; i++) {
        c[i] = 0.f;
    }

    a->data[0][0] = 0.31; a->data[0][1] = 0.14; a->data[0][2] = 0.30; a->data[0][3] = 0.27; a->data[1][0] = 0.26; a->data[1][1] = 0.32; a->data[1][2] = 0.18; a->data[1][3] = 0.24; a->data[2][0] = 0.61; a->data[2][1] = 0.22; a->data[2][2] = 0.20; a->data[2][3] = 0.31; a->data[3][0] = 0.4; a->data[3][1] = 0.34; a->data[3][2] = 0.36; a->data[3][3] = 0.17; 
    b[0] = 1.02;
    b[1] = 1.00;
    b[2] = 1.34;
    b[3] = 1.27;

    printf("\nsoal 1\n");
    m_print(a);
    printf("gauss naive: ");
    gaussnaive(a, b, c);
    for (i = 0; i < size; i++) {
        printf("%.5f ", c[i]);
        c[i] = 0.f;
    }
    printf("\n");
    printf("gauss jordan: ");
    gaussjordan(a, b, c);
    for (i = 0; i < size; i++) {
        printf("%.5f ", c[i]);
        c[i] = 0.f;
    }
    printf("\n");
    printf("LU decomposition: ");
    ludecomp(a, b, c);
    for (i = 0; i < size; i++) {
        printf("%.5f ", c[i]);
        c[i] = 0.f;
    }
    printf("\n");

    m_del(a);
    free(b);
    free(c);

    /***** SOAL 2 *****/
    size = 9;
    a = m_init(size, size);
    b = (double *) malloc(size * sizeof(double));
    c = (double *) malloc(size * sizeof(double));

    for (i = 0; i < size; i++) {
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
    b[0] = 0;
    b[1] = -1000;
    b[2] = 0;
    b[3] = 0;
    b[4] = 500;
    b[5] = 0;
    b[6] = 0;
    b[7] = -500;
    b[8] = 0;

    printf("\nsoal 2\n");
    m_print(a);
    printf("gauss naive: ");
    gaussnaive(a, b, c);
    for (i = 0; i < size; i++) {
        printf("%.5f ", c[i]);
        c[i] = 0.f;
    }
    printf("\n");
    printf("gauss jordan: ");
    gaussjordan(a, b, c);
    for (i = 0; i < size; i++) {
        printf("%.5f ", c[i]);
        c[i] = 0.f;
    }
    printf("\n");
    printf("LU decomposition: ");
    ludecomp(a, b, c);
    for (i = 0; i < size; i++) {
        printf("%.5f ", c[i]);
        c[i] = 0.f;
    }
    printf("\n");

    m_del(a);
    free(b);
    free(c);

    /**SOAL 3**/
    size = 6;
    a = m_init(size, size);
    b = (double *) malloc(size * sizeof(double));
    c = (double *) malloc(size * sizeof(double));

    for (i = 0; i < size; i++) {
        b[i] = 1.f;
        c[i] = 0.f;
    }

    printf("\n soal 3\n");
    m_hilbert(a);
    m_print(a);
    ludecomp(a, b, c);
    //m_print(a);
    printf("LU decomposition: ");
    ludecomp(a, b, c);
    for (i = 0; i < size; i++) {
        printf("%.5f ", c[i]);
        c[i] = 0.f;
    }
    printf("\n");

    m_del(a);
    free(b);
    free(c);

    /** SOAL 4 (TERAPAN) **/
    size = 6;
    a = m_init(size, size);
    b = (double *) malloc(size * sizeof(double));
    c = (double *) malloc(size * sizeof(double));

    for (i = 0; i < size; i++) {
        c[i] = 0.f;
    }

    a->data[0][0] = 0.866; a->data[0][1] = 0; a->data[0][2] = -0.5; a->data[0][3] = 0; a->data[0][4] = 0; a->data[0][5] = 0; a->data[1][0] = 0.5; a->data[1][1] = 0; a->data[1][2] = 0.866; a->data[1][3] = 0; a->data[1][4] = 0; a->data[1][5] = 0; a->data[2][0] = -0.866; a->data[2][1] = -1; a->data[2][2] = 0; a->data[2][3] = -1; a->data[2][4] = 0; a->data[2][5] = 0; a->data[3][0] = -0.5; a->data[3][1] = 0; a->data[3][2] = 0; a->data[3][3] = 0; a->data[3][4] = -1; a->data[3][5] = 0; a->data[4][0] = 0; a->data[4][1] = 1; a->data[4][2] = 0.5; a->data[4][3] = 0; a->data[4][4] = 0; a->data[4][5] = 0; a->data[5][0] = 0; a->data[5][1] = 0; a->data[5][2] = -0.866; a->data[5][3] = 0; a->data[5][4] = 0; a->data[5][5] = -1; 
    b[0] = 0;
    b[1] = -1000;
    b[2] = 0;
    b[3] = 0;
    b[4] = 0;
    b[5] = 0;

    printf("\nsoal terapan\n");
    m_print(a);
    printf("LU decomposition: ");
    ludecomp(a, b, c);
    for (i = 0; i < size; i++) {
        printf("%.5f ", c[i]);
        c[i] = 0.f;
    }
    printf("\n");

    m_del(a);
    free(b);
    free(c);

    return 0;
}
