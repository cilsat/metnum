#include <stdlib.h>
#include <stdio.h>

#define PREC "%.3f\t"
#define MAX_VAL 100

typedef struct {
    long rows;
    long cols;
    long mem_size;
    double *data;
} matrix;

// Initializes a matrix of specified rows and columns size. The data field
// is a 2-d contiguous array in memory.
matrix *m_init(long rows, long cols) {
    matrix *m = (matrix *)malloc(sizeof(matrix));

    m->rows = rows;
    m->cols = cols;
    m->mem_size = rows*cols*sizeof(double);

    double *m_temp = (double *)malloc(rows*cols*sizeof(double));
    for (long i = 0; i < rows; i++) {
        for (long j = 0; j < cols; j++) {
            m_temp[i*cols + j] = 0.f;
        }
    }
    m->data = m_temp;

    return m;
}

// Destructor for square matrix
void m_del(matrix *m) {
    free(m->data);
    free(m);
}

void m_print(matrix *m) {
    for (long i = 0; i < m->rows; i++) {
        for (long j = 0; j < m->cols; j++) {
            printf(PREC, m->data[i*m->cols + j]);
        }
        printf("\n");
    }
    printf("\n");
}

matrix *m_slice(matrix *m, long r_start, long r_stop, long c_start, long c_stop) {
    matrix *temp = m_init(r_stop-r_start, c_stop-c_start);
    for (long i = r_start; i < r_stop; i++) {
        for (long j = c_start; j < c_stop; j++) {
            temp->data[(i-r_start)*c_stop + (j-c_start)] = m->data[i*m->cols + j];
        }
    }

    return temp;
}

// Calculates the dot product of A and B and stores results in C.
__global__ void m_dot(double *c, double *a, double *b, long a_rows, long a_cols, long b_rows, long b_cols) {
    long i, j, k;
    double sum = 0;

    if (a_cols != b_rows) {
        printf("these matrices cannot be multiplied!\n");
    }

    for (i = 0; i < a_rows; i++) {
        for (j = 0; j < b_cols; j++) {
            sum = 0;
            for (k = 0; k < a_cols; k++) {
                sum += a[i*a_cols + k]*b[k*b_rows + j];
            }
            c[i*b_cols + j] = sum;
        }
    }
}

// Randomizes the data field of a given matrix
void m_rand(matrix *m) {
    long i, j;
    srand(time(NULL));

    for (i = 0; i < m->rows; i++) {
        for (j = 0; j < m->cols; j++) {
            m->data[i*m->cols + j] = MAX_VAL*((double) rand() / (double) RAND_MAX);
        }
    }
}

__global__ void m_hilbert(double *m, long size) {
    long n = size;
    long i, j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            m[i*n + j] = 1.f/(i + j + 1);
        }
    }
}

