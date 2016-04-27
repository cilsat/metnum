#include "num.h"
#include <string.h>
#include <math.h>

#define INF 99999.9
#define MAX_ITER 1000
#define MAX_ERR 0.0001

void gaussjordan(matrix *a, double *b, double* results) {
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

    // solve
    for (i = 0; i < m->rows; i++) {
        results[i] = m->data[i][m->cols-1];
    }
    m_del(m);
}

void gaussnaive(matrix *a, double *b, double *results) {
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
            sum -= m->data[i][j]*results[j];
        }
        results[i] = sum/m->data[i][i];
    }
    // free matrix
    m_del(m);
}

void ludecomp(matrix *a, double *b, double *results) {
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
            double n_pivot = u->data[j][i]/u->data[i][i];
            l->data[j][i] = n_pivot;
            for (k = 0; k < u->cols; k++) {
                u->data[j][k] -= n_pivot * u->data[i][k];
            }
        }
    }

    // forward substitution
    double *d = (double *) malloc(size*sizeof(double));
    for (i = 0; i < size; i++) {
        d[i] = 0.f;
    }

    for (i = 0; i < size; i++) {
        double sum = b[i];
        for (j = 0; j <= i; j++) {
            sum -= l->data[i][j] * d[j];
        }
        d[i] = sum/l->data[i][i];
    }

    m_print(u);
    // back substitution
    for (i = size-1; i >=0; i--) {
        double sum = d[i];
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
    double gx = (70 + 1.463/pow(x,2))*(x-0.0394) - 0.08314*215;
    return gx;
}

double bisection(nirlanjar f, double xmin, double xmax) {
    double xmid, xmid_prev, fxmin, fxmax, test, root;
    double ea = INF;
    double es = MAX_ERR;
    long iter = 0;
    
    fxmin = f(xmin);
    fxmax = f(xmax);
    xmid = xmin;

    if (fxmin*fxmax > 0) {
        printf("rentang tidak valid\n");
        return 0;
    }
    
    while (ea > es && iter < MAX_ITER) {
        xmid_prev = xmid;
        xmid = 0.5*(xmin + xmax);
        ea = fabs((xmid - xmid_prev)/xmid)*100;
        iter++;
        test = f(xmid)*f(xmin);
        if (test < 0) 
            xmax = xmid;
        else
            xmin = xmid;
	}
    root = xmid;
    return root;
}

double falsepos(nirlanjar f, double xmin, double xmax) {
    double xl, xu, xr, fxmin, fxmax;
    xl = xmin;
    xu = xmax;
    long iter = 0;

    fxmin = f(xl);
    fxmax = f(xu);

    if (fxmin*fxmax > 0) {
        printf("rentang tidak valid\n");
        return 0;
    }

    xr = xu - (fxmax*(xl - xu)/(fxmin - fxmax));
    while (fabs(xl - xu) > MAX_ERR && iter < MAX_ITER) {
        iter++;
        if (fxmin*f(xr) < 0)
            xu = xr;
        else
            xl = xr;
        xr = xu - (f(xu)*(xl-xu)/(f(xl)-f(xu)));
    }
    return xr;
}

double secant(nirlanjar f, double x1, double x2) {
    double x3, fx1, fx2, er, eps;
    int iter = 0;

    fx1 = f(x1);
    fx2 = f(x2);
    x3 = x2 - fx2*(x1 - x2)/(fx1 - fx2);

    er = fabs(x3 - x2);
    eps = 2*er/(fabs(x3) + MAX_ERR);
    while ((er > MAX_ERR) | (eps > MAX_ERR) && (iter < MAX_ITER)) {
        iter++;
        x1 = x2;
        x2 = x3;
        fx1 = f(x1);
        fx2 = f(x2);
        x3 = x2 - fx2*(x1 - x2)/(fx1 - fx2);
        er = fabs(x3 - x2);
        eps = 2*er/(fabs(x3) + MAX_ERR);
    }
    return x3;
}

double newton(nirlanjar f, nirlanjar df, double x) {
    double xr, xn, ea;
    ea = INF;
    xr = x;
    int iter = 0;

    xn = xr - f(xr)/df(xr);
    while (ea > MAX_ERR && iter < MAX_ITER) {
        iter++;
        if (xn != 0)
            ea = fabs((xn - xr)/xr)*100;
        xr = xn;
        xn = xr - f(xr)/df(xr);
    }
    return xn;
}
