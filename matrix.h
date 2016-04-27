#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAX_double 100

typedef struct {
    long rows;
    long cols;
    long size;
    double **data;
} matrix;

matrix *m_init(long, long);
void m_rand(matrix *);
void m_del(matrix *);
void m_print(matrix *);
void m_div(matrix *, double);
matrix *m_mul(matrix *, double);
void m_add(matrix *, double);
void m_sub(matrix *, double);
matrix *m_slic(matrix *, long, long, long, long);
void m_dot(matrix *, matrix *, matrix *);
void m_hilbert(matrix *);

#endif
