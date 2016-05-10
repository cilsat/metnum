#include "matrix.cu"
#include <string.h>

#define GRID_SIZE 100
#define BLOCK_SIZE 512
#define MATRIX_SIZE 10

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

    double *hA = (double *)malloc(msize);
    double *hb = (double *)malloc(vsize);
    double *hc = (double *)malloc(vsize);
    double *dA, *db, *dc;

    cudaMalloc((void **) &dA, msize);
    cudaMalloc(&db, vsize);
    cudaMalloc(&dc, vsize);

    cudaMemcpy(dA, hA, msize, cudaMemcpyHostToDevice);
    cudaMemcpy(db, hb, vsize, cudaMemcpyHostToDevice);
    cudaMemcpy(dc, hc, vsize, cudaMemcpyHostToDevice);

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    /*
    h = m_init(n, n);
    cudaMallocPitch(&d, &pitch, n*sizeof(double), n);
    cudaMalloc(&dest, n*n*sizeof(double));

    cu_hilbert<<< GRID_SIZE, BLOCK_SIZE >>>(d, dest, pitch, n, n);
    cudaMemcpy(h->data, dest, n*n*sizeof(double), cudaMemcpyDeviceToHost);
    m_print(h);

    m_del(h);
    cudaFree(d);
    cudaFree(dest);
    */

    cudaFree(dA);
    cudaFree(db);
    cudaFree(dc);
    free(hA);
    free(hb);
    free(hc);
    return 0;
}
