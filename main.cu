#include "matrix.cu"
#include "num.cu"
#include <string.h>

#define GRID_SIZE 9
#define BLOCK_SIZE 9
#define MATRIX_SIZE 3

// enumerate function names
enum functions {
    naive = 0,
    jordan,
    doolittle,
    crout,
    jacobi,
    bisection,
    falsepos,
    secant
};

// matrix and numerical function pointers
typedef void (*m_func) (double *A, double *b, double *c, long n);
typedef void (*n_func) (function f, double x0, double x1);

__global__ void m_calc(m_func f, char *file) {
    // read file
    
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
    char filename[16] = "something";
    m_calc(crout, "something");

    return 0;
}
