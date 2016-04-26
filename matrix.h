#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAX_FLOAT 100

typedef struct {
    long rows;
    long cols;
    long size;
    float **data;
} matrix;

matrix *m_init(long, long);
void m_rand(matrix *);
void m_del(matrix *);
void m_print(matrix *);
void m_div(matrix *, float);
matrix *m_mul(matrix *, float);
void m_add(matrix *, float);
void m_sub(matrix *, float);
matrix *m_slic(matrix *, long, long, long, long);
void m_dot(matrix *, matrix *, matrix *);

#endif
