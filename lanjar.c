#include "num.h"
#include <stdio.h>

int main() {
    long i;

    /***** SOAL 1 *****/
    long size = 4;
    matrix *a = m_init(size, size);
    float *b = (float *) malloc(size * sizeof(float));
    float *c = (float *) malloc(size * sizeof(float));

    for (i = 0; i < size; i++) {
        c[i] = 0.f;
    }

    a->data[0][0] = 0.31;
    a->data[0][1] = 0.14;
    a->data[0][2] = 0.30;
    a->data[0][3] = 0.27;
    a->data[1][0] = 0.26;
    a->data[1][1] = 0.32;
    a->data[1][2] = 0.18;
    a->data[1][3] = 0.24;
    a->data[2][0] = 0.61;
    a->data[2][1] = 0.22;
    a->data[2][2] = 0.20;
    a->data[2][3] = 0.31;
    a->data[3][0] = 0.4;
    a->data[3][1] = 0.34;
    a->data[3][2] = 0.36;
    a->data[3][3] = 0.17;

    b[0] = 1.02;
    b[1] = 1.00;
    b[2] = 1.34;
    b[3] = 1.27;


    /***** SOAL 2 *****/


    return 0;
}
