#include "matrix.h"

#define PREC "%.3f\t"

// Initializes a matrix of specified rows and columns size. The data field
// is a 2-d contiguous array in memory.
matrix *m_init(long rows, long cols) {
    matrix *m = calloc(1, sizeof(matrix));

    m->rows = rows;
    m->cols = cols;
    m->size = rows*cols;

    double **m_temp = (double **)malloc(rows*sizeof(double *));
    double *data = (double *)malloc(rows*cols*sizeof(double));
    for (long i = 0; i < rows; i++) {
        m_temp[i] = &(data[cols*i]);
        for (long j = 0; j < cols; j++) {
            m_temp[i][j] = 0.f;
        }
    }
    m->data = m_temp;

    return m;
}

// Destructor for square matrix
void m_del(matrix *m) {
    free(m->data[0]);
    free(m->data);
    free(m);
}

// Deep copy matrix src to dest. This is needed because our matrix struct
// contains a multi-dimensional array, and hence shallow copy with memcpy
// only copies the pointers but not what they're pointing to.
void m_cpy(matrix *dest, matrix *src) {
    long i, j;

    // check size
    if (dest->size != src->size)
        printf("source and destinations matrix are of different size!\n");

    // copy non-data fields directly
    dest->rows = src->rows;
    dest->cols = src->cols;
    dest->size = src->size;

    // copy data field 
    for (i = 0; i < src->rows; i++) {
        for (j = 0; j < src->cols; j++) {
            dest->data[i][j] = src->data[i][j];
        }
    }
}

void m_print(matrix *m) {
    for (long i = 0; i < m->rows; i++) {
        for (long j = 0; j < m->cols; j++) {
            printf(PREC, m->data[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

matrix *m_slice(matrix *m, long r_start, long r_stop, long c_start, long c_stop) {
    matrix *temp = m_init(r_stop-r_start, c_stop-c_start);
    for (long i = r_start; i < r_stop; i++) {
        for (long j = c_start; j < c_stop; j++) {
            temp->data[i-r_start][j-c_start] = m->data[i][j];
        }
    }

    return temp;
}

// Calculates the dot product of A and B and stores results in C.
void m_dot(matrix *a, matrix *b, matrix *c) {
    long i, j, k;
    double sum = 0;

    if (a->cols != b->rows) {
        printf("these matrices cannot be multiplied!\n");
    } 
    if ((c->rows != a->rows) | (c->cols != b->cols)) {
        printf("size of output matrix invalid!\n");
    }

    for (i = 0; i < c->rows; i++) {
        for (j = 0; j < c->cols; j++) {
            sum = 0;
            for (k = 0; k < a->cols; k++) {
                sum += a->data[i][k]*b->data[k][j];
            }
            c->data[i][j] = sum;
        }
    }
}

// Randomizes the data field of a given matrix
void m_rand(matrix *m) {
    long i, j;
    srand(time(NULL));

    for (i = 0; i < m->rows; i++) {
        for (j = 0; j < m->cols; j++) {
            m->data[i][j] = MAX_VAL*((double) rand() / (double) RAND_MAX);
        }
    }
}

void m_hilbert(matrix *m) {
    long n = m->rows;
    long i, j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            m->data[i][j] = 1.f/(i + j + 1);
        }
    }
}
