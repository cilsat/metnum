#include <string.h>
#include <math.h>

#define INF 99999
#define MAX_ITER 100
#define MAX_ERR 0.0001
#define INIT_JACOBI 0
#define DEBUG 1

__device__ void pivot(double *m, long m_rows, long m_cols, long cur_row) {
    long i = cur_row;
    long j;
    double temp = 0;
    double *d;

    d = (double *) malloc(m_cols*sizeof(double));

    double p = m[i*m_cols + i];
    if (p == 0) {
        // where do we go?
        long new_row = i;
        double max_row = 0;
        for (j = i+1; j < m_rows; j++) {
            temp = fabs(m[j*m_cols + i]);
            if ((i != j) && (temp*m[i*m_cols + j] != 0) && (temp >= max_row)) {
                max_row = temp;
                new_row = j;
            }
        }
        // swap cur_row with new_row
        for (j = 0; j < m_cols; j++) {
            d[j] = m[new_row*m_cols + j];
            m[new_row*m_cols + j] = m[i*m_cols + j];
            m[i*m_cols + j] = d[j];
        }
    }
    free(d);
}

__global__ void gaussjordan(double *a, double *b, double* c, long a_rows) {
    long i, j, k;

    // combine left and right sides
    double *m;
    long m_rows, m_cols;
    m_rows = a_rows;
    m_cols = a_rows + 1;
    m = (double *) malloc(m_rows*m_cols*sizeof(double));

    for (i = 0; i < m_rows; i++) {
        for (j = 0; j < a_rows; j++) {
            m[i*a_rows + j] = a[i*a_rows + j];
        }
        m[i*a_rows + a_rows] = b[i];
    }

    for (i = 0; i < m_rows; i++) {
        // pivot
        pivot(m, m_rows, m_cols, i);
        // normalize row
        double pivot = m[i*m_cols + i];
        for (k = 0; k < m_cols; k++) {
            m[i*m_cols + k] /= pivot;
        }
        // subtract row from all other rows
        for (j = 0; j < m_rows; j++) {
            if (j != i) {
                double n_pivot = m[j*m_cols + i];
                for (k = 0; k < m_cols; k++) {
                    m[j*m_cols + k] -= n_pivot * m[i*m_cols + k];
                }
            }
        }
    }

    // forward substitution
    for (i = 0; i < m_rows; i++) {
        c[i] = m[i*m_cols + m_cols-1];
    }
    free(m);
}

__global__ void gaussnaive(double *a, double *b, double *c, long a_rows) {
    long i, j, k;

    // combine left and right sides
    double *m, *d;
    long m_rows = a_rows;
    long m_cols = a_rows + 1;
    m = (double *) malloc(m_rows*m_cols*sizeof(double));
    d = (double *) malloc(m_rows*sizeof(double));

    for (i = 0; i < m_rows; i++) {
        for (j = 0; j < m_rows; j++) {
            m[i*m_cols + j] = a[i*m_rows + j];
        }
        m[i*m_cols + m_rows] = b[i];
        d[i] = 0.f;
    }

    // forward elimination
    for (i = 0; i < m_rows; i++) {
        // pivot
        pivot(m, m_rows, m_cols, i);
        // subtract row from all subsequent rows
        for (j = i+1; j < m_rows; j++) {
            double n_pivot = m[j*m_cols + i]/m[i*m_cols + i];
            for (k = 0; k < m_cols; k++) {
                m[j*m_cols + k] -= n_pivot * m[i*m_cols + k];
            }
        }
    }

    // back substitution
    for (i = m_rows-1; i >=0; i--) {
        double sum = m[i*m_cols + m_rows];
        for (j = i; j < m_rows; j++) {
            sum -= m[i*m_cols + j]*d[j];
        }
        c[i] = sum/m[i*m_cols + i];
        d[i] = c[i];
    }

    if (DEBUG) {
        for (i = 0; i < m_rows; i++) {
            for (j = 0; j < m_rows; j++) {
                a[i*a_rows + j] = m[i*m_cols + j];
            }
            b[i] = m[i*m_cols + m_rows];
        }
    }

    // free matrix
    free(m);
    free(d);
}

__global__ void doolittle(double *a, double *b, double *c, long a_rows) {
    long i, j, k;
    long n = a_rows;

    // create and initialize lower, upper, and intermediate results  matrixes
    double *l, *u, *d, *temp;
    l = (double *) malloc(n*n*sizeof(double));
    u = (double *) malloc(n*(n+1)*sizeof(double));
    d = (double *) malloc(n*sizeof(double));
    temp = (double *) malloc(n*sizeof(double));

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            u[i*(n+1) + j] = a[i*n + j];
        }
        u[i*(n+1) + n] = b[i];
        l[i*n + i] = 1.f;
        d[i] = 0.f;
        temp[i] = 0.f;
    }

    // decomposition
    for (i = 0; i < n; i++) {
        // pivot
        pivot(u, n, n+1, i);
        // subtract row from all subsequent rows
        for (j = i+1; j < n; j++) {
            double n_pivot = u[j*(n+1) + i]/u[i*(n+1) + i];
            l[j*n + i] = n_pivot;
            for (k = 0; k < n; k++) {
                u[j*(n+1) + k] -= n_pivot * u[i*(n+1) + k];
            }
        }
    }

    // forward substitution
    for (i = 0; i < n; i++) {
        double sum = u[i*(n+1) + n];
        for (j = 0; j <= i; j++) {
            sum -= l[i*n + j] * temp[j];
        }
        d[i] = sum/l[i*n + i];
        temp[i] = d[i];
    }

    for (i = 0; i < n; i++) {
        temp[i] = 0.f;
    }

    // back substitution
    for (i = n-1; i >=0; i--) {
        double sum = d[i];
        for (j = i; j < n; j++) {
            sum -= u[i*(n+1) + j] * temp[j];
        }
        c[i] = sum/u[i*(n+1) + i];
        temp[i] = c[i];
    }

    free(l);
    free(u);
    free(d);
    free(temp);
}

__global__ void crout(double *a, double *b, double *c, long a_rows) {
    long i, j, k;
    long n = a_rows;
    double sum;

    double *l, *u, *d;
    l = (double *) malloc(n*n*sizeof(double));
    u = (double *) malloc(n*(n+1)*sizeof(double));
    d = (double *) malloc(n*sizeof(double));

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            u[i*n + j] = a[i*n + j];
        }
        u[i*n + n] = b[i];
        l[i*n + i] = 1.f;
        d[i] = 0.f;
    }

    // decomposition
    for (i = 0; i < n; i++) {
        // pivot
        pivot(u, n, n+1,  i);
        // normalize
        double f_pivot = u[i*n + i];
        l[i*n + i] = f_pivot;
        for (k = i; k < n; k++) {
            u[i*n + k] /= f_pivot;
        }
        // subtract for all subsequent rows
        for (j = i+1; j < n; j++) {
            double n_pivot = u[j*n + i];
            l[j*n + i] = n_pivot;
            for (k = 0; k < n; k++) {
                u[j*n + k] -= n_pivot * u[i*n + k];
            }
        }
    }

    // forward substitution
    for (i = 0; i < n; i++) {
        sum = u[i*n + n];
        for (j = 0; j <= i; j++) {
            sum -= l[i*n + j] * d[j];
        }
        d[i] = sum/l[i*n + i];
    }

    // back substitution
    for (i = n-1; i >=0; i--) {
        sum = d[i];
        for (j = i; j < n; j++) {
            sum -= u[i*n + j] * c[j];
        }
        c[i] = sum/u[i*n + i];
    }
    free(l);
    free(u);
    free(d);
}

__global__ void jacobi(double *a, double *b, double *c, long a_rows) {
    long i, k;
    long n = a_rows;
    int iter = 0;
    double tot_err = INF;
    double sum;

    double *m, *c_prev;
    m = (double *) malloc(n*(n+1)*sizeof(double));
    c_prev = (double *) malloc(n*sizeof(double));

    for (i = 0; i < n; i++) {
        for (k = 0; k < n; k++) {
            m[i*n + k] = a[i*n + k];
        }
        m[i*n + n] = b[i];
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
                    sum += c[k] * m[i*(n+1) + k];
            }
            c[i]= (m[i*(n+1) + n] - sum)/m[i*(n+1) + i];
            if (c[i] != 0)
                tot_err += fabs((c[i] - c_prev[i])/c[i]);
            c_prev[i] = c[i];
        }
        tot_err /= n;
    }
    free(m);
    free(c_prev);
}

typedef double (*function) (double x);

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
