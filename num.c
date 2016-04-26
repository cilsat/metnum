#include "num.h"
#include <string.h>

void gaussjordan(matrix *a, float *b, float* results) {
    long i, j, k;

    // combine left and right sides
    matrix *m = m_init(a->rows, a->cols+1);
    for (i = 0; i < m->rows; i++) {
        for (j = 0; j < a->cols; j++) {
            m->data[i][j] = a->data[i][j];
        }
        m->data[i][m->cols-1] = b[i];
    }

    for (i = 0; i < m->rows; i++) {
        // normalize row
        float pivot = m->data[i][i];
        for (k = 0; k < m->cols; k++) {
            m->data[i][k] /= pivot;
        }

        // subtract row from all other rows
        for (j = 0; j < m->rows; j++) {
            if (j != i) {
                float n_pivot = m->data[j][i];
                for (k = 0; k < m->cols; k++) {
                    m->data[j][k] -= n_pivot * m->data[i][k];
                }
            }
        }
    }

    // solve
    for (i = 0; i < m->rows; i++) {
        results[i] = m->data[i][m->cols-1];
    }
    m_del(m);
}

void gaussnaive(matrix *a, float *b, float *results) {
    long i, j, k;

    // combine left and right sides
    matrix *m = m_init(a->rows, a->cols+1);
    for (i = 0; i < m->rows; i++) {
        for (j = 0; j < a->cols; j++) {
            m->data[i][j] = a->data[i][j];
        }
        m->data[i][m->cols-1] = b[i];
    }

    // forward elimination
    for (i = 0; i < m->rows; i++) {
        // subtract row from all subsequent rows
        for (j = i+1; j < m->rows; j++) {
            float n_pivot = m->data[j][i]/m->data[i][i];
            for (k = 0; k < m->cols; k++) {
                m->data[j][k] -= n_pivot * m->data[i][k];
            }
        }
    }

    // back substitution
    for (i = m->rows-1; i >=0; i--) {
        float sum = m->data[i][m->cols-1];
        for (j = i; j < m->cols-1; j++) {
            sum -= m->data[i][j]*results[j];
        }
        results[i] = sum/m->data[i][i];
    }
    // free matrix
    m_del(m);
}

void ludecomp(matrix *a, float *b, float *results) {
    long i, j, k;
    long size = a->rows;

    matrix *u = m_init(a->rows, a->cols);
    memcpy(u, a, sizeof(*a));
    matrix *l = m_init(a->rows, a->cols);
    for (i = 0; i < size; i++) {
        for (j = 0; j < l->cols; j++) {
            if (i == j)
                l->data[i][j] = 1.f;
            else
                l->data[i][j] = 0.f;
        }
    }

    // decomposition
    for (i = 0; i < u->rows; i++) {
        // subtract row from all subsequent rows
        for (j = i+1; j < u->rows; j++) {
            float n_pivot = u->data[j][i]/u->data[i][i];
            l->data[j][i] = n_pivot;
            for (k = 0; k < u->cols; k++) {
                u->data[j][k] -= n_pivot * u->data[i][k];
            }
        }
    }

    // forward substitution
    float *d = (float *) malloc(size*sizeof(float));
    for (i = 0; i < size; i++) {
        d[i] = 0.f;
    }

    for (i = 0; i < size; i++) {
        float sum = b[i];
        for (j = 0; j <= i; j++) {
            sum -= l->data[i][j] * d[j];
        }
        d[i] = sum/l->data[i][i];
    }

    m_print(u);
    // back substitution
    for (i = size-1; i >=0; i--) {
        float sum = d[i];
        for (j = i; j < size; j++) {
            sum -= u->data[i][j] * results[j];
        }
        results[i] = sum/u->data[i][i];
    }
    for (i = 0; i < size; i++) {
        printf("%.6f ", results[i]);
    }

    free(d);
    m_del(l);
}

