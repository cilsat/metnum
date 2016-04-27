#include "matrix.h"

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

void m_rand(matrix *m) {
    long i, j;
    srand(time(NULL));

    for (i = 0; i < m->rows; i++) {
        for (j = 0; j < m->cols; j++) {
            m->data[i][j] = MAX_double*((double) rand() / (double) RAND_MAX);
        }
    }
}

// Destructor for square matrix
void m_del(matrix *m) {
    free(m->data[0]);
    free(m->data);
    free(m);
}

void m_print(matrix *m) {
    for (long i = 0; i < m->rows; i++) {
        for (long j = 0; j < m->cols; j++) {
            printf("%f ", m->data[i][j]);
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

void m_dot(matrix *a, matrix *b, matrix *c) {
    long i, j;

    if (a->cols != b->rows) {
        printf("these matrices cannot be multiplied!\n");
    } 
    for (i = 0; i < c->rows; i++) {
        for (j = 0; j < c->cols; j++) {
            c->data[i][j] += a->data[i][j]*b->data[j][i];
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
