#include "num.h"

#define MAX_double 100
#define DEBUG 2
#define PREC "%.3f "

void reset(double *arr, long n) {
    for (long i = 0; i < n; i++) {
        printf(PREC, arr[i]);
        arr[i] = 0.f;
    }
}

int main(int argc, char *argv[]) {
    long i, n;
    matrix *a;
    double *b, *c;

    /*
    n = 3;
    a = m_init(n, n);
    b = (double *) malloc(n*sizeof(double));
    c = (double *) malloc(n*sizeof(double));
    for (i = 0; i < n; i++) {
        c[i] = 0.f;
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

    printf("initial matrix:\n");
    m_print(a);

    printf("\ncrout\n");
    crout(a, b, c);
    reset(c, n);

    printf("\ndoolittle\n");
    doolittle(a, b, c);
    reset(c, n);

    m_del(a);
    free(b);
    free(c);
    */

    n = strtol(argv[1], NULL, 10);
    a = m_init(n, n);
    b = (double *)malloc(n*sizeof(double));
    c = (double *)malloc(n*sizeof(double));
    for (i = 0; i < n; i++) {
        b[i] = 1.f;
        c[i] = 0.f;
    }

    printf("hilbert matrix size %ld:\n", n);
    m_hilbert(a);
    m_print(a);

    printf("\ncrout\n");
    crout(a, b, c);
    reset(c, n);

    printf("\ndoolittle\n");
    doolittle(a, b, c);
    reset(c, n);

    m_del(a);
    free(b);
    free(c);

    return 0;
}
