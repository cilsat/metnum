#include "num.h"
#include <omp.h>

#define MAX_FLOAT 100
#define DEBUG 2

int main(int argc, char *argv[]) {
    //long size = (long) strtol(argv[1], NULL, 10);
    long i;
    double dstart, dstop;

    /*
    matrix *a = m_init(size, size);
    float *b = (float *) malloc(size*sizeof(float));
    float *sol = (float *) malloc(a->rows * sizeof(float));

    m_rand(a);
    srand(time(NULL));
    for (i = 0; i < size; i++) {
        b[i] = MAX_FLOAT*((float) rand() / (float) RAND_MAX);
        sol[i] = 0.f;
    }
    */

    matrix *a = m_init(3, 3);
    float *b = (float *) malloc(3*sizeof(float));
    float *sol = (float *) malloc(a->rows * sizeof(float));
    for (i = 0; i < 3; i++) {
        sol[i] = 0.f;
    }

    a->data[0][0] = 3.0;
    a->data[0][1] = -0.1;
    a->data[0][2] = -0.2;
    a->data[1][0] = 0.1;
    a->data[1][1] = 7.0;
    a->data[1][2] = -0.3;
    a->data[2][0] = 0.3;
    a->data[2][1] = -0.2;
    a->data[2][2] = 10.0;

    b[0] = 7.85;
    b[1] = -19.3;
    b[2] = 71.4;

    dstart = omp_get_wtime();
    gaussnaive(a, b, sol);
    dstop = omp_get_wtime();

    if (DEBUG == 2) {
        printf("%.5f\n", dstop-dstart);
    }

    if (DEBUG == 1) {
        m_print(a);
        for (i = 0; i < a->rows; i++) {
            printf("%f ", sol[i]);
        }
    }
    for (i = 0; i < 3; i++) {
        sol[i] = 0.f;
    }

    dstart = omp_get_wtime();
    gaussjordan(a, b, sol);
    dstop = omp_get_wtime();

    if (DEBUG == 2) {
        printf("%.5f\n", dstop-dstart);
    }
    for (i = 0; i < 3; i++) {
        sol[i] = 0.f;
    }

    dstart = omp_get_wtime();
    ludecomp(a, b, sol);
    dstop = omp_get_wtime();

    if (DEBUG == 2) {
        printf("%.5f\n", dstop-dstart);
    }

    m_del(a);
    free(b);
    free(sol);
    return 0;
}
