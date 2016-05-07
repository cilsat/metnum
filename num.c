#include "num.h"
#include <string.h>
#include <math.h>

#define INF 99999
#define MAX_ITER 100
#define MAX_ERR 0.0001
#define INIT_JACOBI 0
#define PREC "%.5f "
#define DEBUG 0

void pivot(matrix *m, long cur_row) {
    long i = cur_row;
    long j;
    double temp = 0;
    double *d = (double *) malloc(m->cols * sizeof(double));

    double p = m->data[i][i];
    if (p == 0) {
        // where do we go?
        long new_row = -1;
        double max_row = 0;
        for (j = i+1; j < m->rows; j++) {
            temp = fabs(m->data[j][i]);
            if ((i != j) && (temp*m->data[i][j] != 0) && (temp >= max_row)) {
                max_row = temp;
                new_row = j;
            }
        }
        // swap cur_row with new_row
        if (new_row == -1) {
            printf("matriks tidak memiliki solusi\n");
        } else {
            //printf("nr %ld\n", new_row);
            for (j = 0; j < m->cols; j++) {
                d[j] = m->data[new_row][j];
                m->data[new_row][j] = m->data[i][j];
                m->data[i][j] = d[j];
            }
        }
    }
    free(d);
}

void gaussjordan(matrix *a, double *b, double* c) {
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
        // pivot
        pivot(m, i);
        // normalize row
        double pivot = m->data[i][i];
        for (k = 0; k < m->cols; k++) {
            m->data[i][k] /= pivot;
        }
        // subtract row from all other rows
        for (j = 0; j < m->rows; j++) {
            if (j != i) {
                double n_pivot = m->data[j][i];
                for (k = 0; k < m->cols; k++) {
                    m->data[j][k] -= n_pivot * m->data[i][k];
                }
            }
        }
    }

    // forward substitution
    for (i = 0; i < m->rows; i++) {
        c[i] = m->data[i][m->cols-1];
    }
    m_del(m);
}

void gaussnaive(matrix *a, double *b, double *c) {
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
        // pivot
        pivot(m, i);
        // subtract row from all subsequent rows
        for (j = i+1; j < m->rows; j++) {
            double n_pivot = m->data[j][i]/m->data[i][i];
            for (k = 0; k < m->cols; k++) {
                m->data[j][k] -= n_pivot * m->data[i][k];
            }
        }
    }

    // back substitution
    for (i = m->rows-1; i >=0; i--) {
        double sum = m->data[i][m->cols-1];
        for (j = i; j < m->cols-1; j++) {
            sum -= m->data[i][j]*c[j];
        }
        c[i] = sum/m->data[i][i];
    }
    // free matrix
    m_del(m);
}

void doolittle(matrix *a, double *b, double *c) {
    long i, j, k;
    long n = a->rows;

    // create and initialize lower, upper, and intermediate results  matrixes
    matrix *l = m_init(n, n);
    matrix *u = m_init(n, n+1);
    double *d = (double *)malloc(n*sizeof(double));
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            u->data[i][j] = a->data[i][j];
        }
        u->data[i][n] = b[i];
        l->data[i][i] = 1.f;
        d[i] = 0.f;
    }

    // decomposition
    for (i = 0; i < n; i++) {
        // pivot
        pivot(u, i);
        // subtract row from all subsequent rows
        for (j = i+1; j < n; j++) {
            double n_pivot = u->data[j][i]/u->data[i][i];
            l->data[j][i] = n_pivot;
            for (k = 0; k < n; k++) {
                u->data[j][k] -= n_pivot * u->data[i][k];
            }
        }
    }

    // forward substitution
    for (i = 0; i < n; i++) {
        double sum = u->data[i][n];
        for (j = 0; j <= i; j++) {
            sum -= l->data[i][j] * d[j];
        }
        d[i] = sum/l->data[i][i];
    }

    // back substitution
    for (i = n-1; i >=0; i--) {
        double sum = d[i];
        for (j = i; j < n; j++) {
            sum -= u->data[i][j] * c[j];
        }
        c[i] = sum/u->data[i][i];
    }
    m_del(l);
    m_del(u);
    free(d);
}

void crout(matrix *a, double *b, double *c) {
    long i, j, k;
    long n = a->rows;
    double sum;

    matrix *l = m_init(n, n);
    matrix *u = m_init(n, n+1);

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            u->data[i][j] = a->data[i][j];
        }
        u->data[i][n] = b[i];
    }

    for (i = 0; i < n; i++) {
        l->data[i][i] = 1.f;
    }

    // decomposition
    for (i = 0; i < n; i++) {
        // pivot
        pivot(u, i);
        // normalize
        double f_pivot = u->data[i][i];
        l->data[i][i] = f_pivot;
        for (k = i; k < n; k++) {
            u->data[i][k] /= f_pivot;
        }
        // subtract for all subsequent rows
        for (j = i+1; j < n; j++) {
            double n_pivot = u->data[j][i];
            l->data[j][i] = n_pivot;
            for (k = 0; k < n; k++) {
                u->data[j][k] -= n_pivot * u->data[i][k];
            }
        }
    }

    // forward substitution
    double *d = (double *) malloc(n*sizeof(double));
    for (i = 0; i < n; i++) {
        d[i] = 0.f;
    }

    for (i = 0; i < n; i++) {
        sum = u->data[i][n];
        for (j = 0; j <= i; j++) {
            sum -= l->data[i][j] * d[j];
        }
        d[i] = sum/l->data[i][i];
    }

    // back substitution
    for (i = n-1; i >=0; i--) {
        sum = d[i];
        for (j = i; j < n; j++) {
            sum -= u->data[i][j] * c[j];
        }
        c[i] = sum/u->data[i][i];
    }
    m_del(l);
    m_del(u);
    free(d);
}

void jacobi(matrix *a, double *b, double *c) {
    long i, k;
    long n = a->rows;
    int iter = 0;
    double tot_err = INF;
    double sum;

    matrix *m = m_init(n, n+1);
    for (i = 0; i < n; i++) {
        for (k = 0; k < n; k++) {
            m->data[i][k] = a->data[i][k];
        }
        m->data[i][n] = b[i];
    }

    double *c_prev = (double *)malloc(n*sizeof(double));
    for (i = 0; i < n; i++) {
        c[i] = INIT_JACOBI;
        c_prev[i] = INIT_JACOBI;
    }

    while (tot_err > MAX_ERR && iter < MAX_ITER) {
        iter++;
        tot_err = 0;
        for (i = 0; i < n; i++) {
            sum = 0;
            for (k = 0; k < n; k++) {
                if (k != i)
                    sum += c[k] * m->data[i][k];
            }
            c[i]= (m->data[i][n] - sum)/m->data[i][i];
            if (c[i] != 0)
                tot_err += fabs((c[i] - c_prev[i])/c[i]);
            c_prev[i] = c[i];
        }
        tot_err /= n;
    }

    m_del(m);
    free(c_prev);
}

double f1a(double x) {
    double fx = (x*pow(2.1 - 0.5*x, 0.5)/((1-x)*pow(1.1-0.5*x, 0.5)))-3.69;
    return fx;
}

double f1b(double x) {
    double fx = tan(x) - x + 1;
    return fx;
}

double f1c(double x) {
    double fx = 0.5*exp(x/3) - sin(x);
    return fx;
}

double f2(double x) {
    double fx = 80*exp(-2.f*x) + 20*exp(-0.1*x);
    return fx;
}

double df2(double x) {
    double dfx;
    dfx = 160*exp(-2*x) + 2*exp(-0.1*x);
    return(dfx);
}

double f3(double x) {
    double dx = (70 + 1.463/pow(x,2))*(x-0.0394) - 0.08314*215;
    return dx;
}

double df3(double x) {
    double dfx = 70 - (1.463/pow(x,2)) + (0.115284/pow(x,3));
    return dfx;
}

double bisection(function f, double xmin, double xmax) {
    double xmid, xmid_prev, test;
    double ea = INF;
    long iter = 0;
    xmid = xmin;

    if (f(xmin)*f(xmax) > 0) {
        printf("rentang tidak valid ");
        return 0;
    }

    while (ea > MAX_ERR && iter < MAX_ITER) {
        iter++;
        xmid_prev = xmid;
        xmid = 0.5*(xmin + xmax);
        ea = fabs((xmid - xmid_prev)/xmid);
        test = f(xmid)*f(xmin);
        if (test < 0)
            xmax = xmid;
        else
            xmin = xmid;
    }
    return xmid;
}

double falsepos(function f, double xmin, double xmax) {
    double xl, xu, xr;
    xl = xmin;
    xu = xmax;
    long iter = 0;

    if (f(xl)*f(xu) > 0) {
        printf("rentang tidak valid ");
        return 0;
    }

    xr = xu - (f(xu)*(xl - xu)/(f(xl) - f(xu)));
    while (fabs(xl - xu) > MAX_ERR && iter < MAX_ITER) {
        iter++;
        if (f(xl)*f(xr) < 0)
            xu = xr;
        else
            xl = xr;
        xr = xu - (f(xu)*(xl-xu)/(f(xl)-f(xu)));
    }
    return xr;
}

double secant(function f, double x1, double x2) {
    double x3, er, eps;
    int iter = 0;

    x3 = x2 - f(x2)*(x1 - x2)/(f(x1) - f(x2));
    er = fabs(x3 - x2);
    eps = 2*er/(fabs(x3) + MAX_ERR);

    while (((er > MAX_ERR) | (eps > MAX_ERR)) && (iter < MAX_ITER)) {
        iter++;
        x1 = x2;
        x2 = x3;
        x3 = x2 - f(x2)*(x1 - x2)/(f(x1) - f(x2));
        er = fabs(x3 - x2);
        eps = 2*er/(fabs(x3) + MAX_ERR);
    }
    return x3;
}

double newton(function f, function df, double x) {
    double xr, xn, ea;
    ea = INF;
    xr = x;
    int iter = 0;

    xn = xr - f(xr)/df(xr);
    while (ea > MAX_ERR && iter < MAX_ITER) {
        iter++;
        if (xn != 0)
            ea = fabs((xn - xr)/xr);
        xr = xn;
        xn = xr - f(xr)/df(xr);
    }
    return xn;
}
