#include "num.cu"
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#define GRID_SIZE 9
#define BLOCK_SIZE 9
#define MATRIX_SIZE 3
#define METHOD SECANT
#define PI 3.14159
#define X0 0.001
#define X1 2.999*PI
#define FX 0

// matrix and numerical function pointers
typedef void (*m_func) (double *A, double *b, double *c, long n);
typedef void (*n_func) (int m, double x0, double x1);

void m_calc(enum method m, char *file) {
    // read file

    if ((m == 0) | (m == 1) | (m == 2) | (m == 3) | (m == 4)) {
        // function pointer to m_func
        m_func f;

        // initialize matrix size
        long n = MATRIX_SIZE;

        // initialize variables
        long msize = n*n*sizeof(double);
        long vsize = n*sizeof(double);

        double *hA, *hb, *hc, *dA, *db, *dc;

        hA = (double *)malloc(msize);
        hb = (double *)malloc(vsize);
        hc = (double *)malloc(vsize);

        cudaMalloc((void **) &dA, msize);
        cudaMalloc(&db, vsize);
        cudaMalloc(&dc, vsize);

        for (long i = 0; i < n; i++) {
            hc[i] = 0.f;
        }

        hA[0] = 3.0; hA[1] = -0.1; hA[2] = -0.2; hA[3] = 0.1; hA[4] = 7.0; hA[5] = -0.3; hA[6] = 0.3; hA[7] = -0.2; hA[8] = 10.0;

        hb[0] = 7.85; hb[1] = -19.3; hb[2] = 71.4;

        cudaMemcpy(dA, hA, msize, cudaMemcpyHostToDevice);
        cudaMemcpy(db, hb, vsize, cudaMemcpyHostToDevice);
        cudaMemcpy(dc, hc, vsize, cudaMemcpyHostToDevice);

        cudaEvent_t start, stop;
        cudaEventCreate(&start);
        cudaEventCreate(&stop);

        switch(m) {
        case NAIVE:
            f = gaussnaive;
            break;
        case JORDAN:
            f = gaussjordan;
            break;
        case DOOLITTLE:
            f = doolittle;
            break;
        case CROUT:
            f = crout;
            break;
        case JACOBI:
            f = jacobi;
            break;
        }

        f<<< GRID_SIZE, n >>>(dA, db, dc, n);

        cudaMemcpy(hA, dA, msize, cudaMemcpyDeviceToHost);
        cudaMemcpy(hb, db, vsize, cudaMemcpyDeviceToHost);
        cudaMemcpy(hc, dc, vsize, cudaMemcpyDeviceToHost);

        for (long i = 0; i < n; i++) {
            for (long j = 0; j < n; j++) {
                printf("%.3f ", hA[i*n + j]);
            }
            printf("\t%.3f\t%.3f", hb[i], hc[i]);
            printf("\n");
        }

        cudaFree(dA);
        cudaFree(db);
        cudaFree(dc);
        free(hA);
        free(hb);
        free(hc);
    }
    else if ((m == 5) | (m == 6) | (m == 7)) {
        // initialize variables
        n_func f;
        double x0, x1;
        typeof(y) hy;
        x0 = X0;
        x1 = X1;

        switch(m) {
        case BISECT:
            f = bisection;
            printf("Bisection\n");
            break;
        case SECANT:
            f = secant;
            printf("Secant\n");
            break;
        case FALSI:
            f = falsepos;
            printf("False Position\n");
            break;
        }

        for (int i = 0; i < N; i++) {
            printf("x0=%.3f x1=%.3f y=", x0, x1);
            f<<< GRID_SIZE, 1 >>>(i, x0, x1);
            cudaMemcpyFromSymbol(&hy, y, sizeof(hy), 0, cudaMemcpyDeviceToHost);
            printf("%.5f\n", hy);
        }
    }
    else if (m == 8) { // Newton Rhapson
    }

    cudaDeviceReset();
}

int main(int argc, char **argv) {
    // Initialize CUDA device
    printf("[MetNum CUDA Routines]\n");
    cudaDeviceProp deviceProp;
    cudaSetDevice(0);
    cudaGetDeviceProperties(&deviceProp, 0);
    printf("GPU Device %d: \"%s\" with compute capability %d.%d\n\n", 0, deviceProp.name, deviceProp.major, deviceProp.minor);

    // parse command line arguments

    // call relevant function
    enum method m = METHOD;
    char filename[16] = "something";
    m_calc(m, filename);

    return 0;
}
