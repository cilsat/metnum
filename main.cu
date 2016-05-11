#include "matrix.cu"
#include "num.cu"
#include <string.h>

#define GRID_SIZE 9
#define BLOCK_SIZE 9
#define MATRIX_SIZE 3

int main(int argc, char **argv) {
    // Initialize CUDA device
    printf("[MetNum CUDA Routines]\n");
    cudaDeviceProp deviceProp;
    cudaSetDevice(0);
    cudaGetDeviceProperties(&deviceProp, 0);
    printf("GPU Device %d: \"%s\" with compute capability %d.%d\n\n", 0, deviceProp.name, deviceProp.major, deviceProp.minor);

    // Initialize blocks and dimensions
    long n = MATRIX_SIZE;
    long msize = n*n*sizeof(double);
    long vsize = n*sizeof(double);

    double *hA, *hb, *hc, *dA, *db, *dc;

    hA = (double *)malloc(msize);
    hb = (double *)malloc(vsize);
    hc = (double *)malloc(vsize);

    cudaMalloc((void **) &dA, msize);
    cudaMalloc(&db, vsize);
    cudaMalloc(&dc, vsize);

    hA[0] = 3.0;
    hA[1] = -0.1;
    hA[2] = -0.2;
    hA[3] = 0.1;
    hA[4] = 7.0;
    hA[5] = -0.3;
    hA[6] = 0.3;
    hA[7] = -0.2;
    hA[8] = 10.0;

    hb[0] = 7.85;
    hb[1] = -19.3;
    hb[2] = 71.4;

    for (long i = 0; i < n; i++) {
        hc[i] = 0.f;
    }

    cudaMemcpy(dA, hA, msize, cudaMemcpyHostToDevice);
    cudaMemcpy(db, hb, vsize, cudaMemcpyHostToDevice);
    cudaMemcpy(dc, hc, vsize, cudaMemcpyHostToDevice);

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    doolittle<<< GRID_SIZE, n >>>(dA, db, dc, n);

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

    return 0;

    /*
    h = m_init(n, n);
    cudaMallocPitch(&d, &pitch, n*sizeof(double), n);
    cudaMalloc(&dest, n*n*sizeof(double));

    cu_hilbert<<< GRID_SIZE, BLOCK_SIZE >>>(d, dest, pitch, n, n);
    cudaMemcpy(h, dest, n*n*sizeof(double), cudaMemcpyDeviceToHost);
    m_print(h);

    m_del(h);
    cudaFree(d);
    cudaFree(dest);
    */
}
