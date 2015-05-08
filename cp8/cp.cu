#include "cp.h"
#include <cuda_runtime.h>
#include <cmath>
#include <stdio.h>
#include <iostream>

#define CHECK_CUDA_ERROR(call) do { \
  cudaError_t result_ = (call); \
  if (result_ != cudaSuccess) { \
  fprintf(stderr, #call " failed: %s\n", \
	  cudaGetErrorString(result_)); \
  exit(1); \
  } \
} while(0)

double get_mean(const float* row, int nx) {
  double sum = 0.0, mean;
  for (int x=0; x<nx; x++) {
    sum += (double) row[x];
  }
  mean = sum / (double) nx;
  return mean;
}

double get_root_square_sum(double* row, int nx) {
  double square_sum = 0, root_square_sum;
  for (int x=0; x<nx; x++) {
    square_sum += pow((double) row[x], 2.0);
  }
  root_square_sum = std::sqrt(square_sum);
  return root_square_sum;
}

__global__ void matrix_multiply(int nx, int ny, double* X, float* result) {
  double r = 0;

  int x = threadIdx.x + blockIdx.x * blockDim.x;
  int y = threadIdx.y + blockIdx.y * blockDim.y;

  if (x >= ny || y >= ny) return;

  for (int i=0; i<nx; i++) {
    r += X[x*nx + i] * X[y*nx + i];
  }

  result[y*ny + x] = r;
}

void correlate(int ny, int nx, const float* data, float* result) {
  double row_mean, row_rss;
  double* X = (double*) malloc(sizeof(double) * nx * ny);
  double* dev_X;
  float* dev_result;

  for (int y=0; y<ny; y++) {
    row_mean = get_mean(&data[y*nx], nx);

    for (int x=0; x<nx; x++) {
      X[y*nx + x] = ((double) data[y*nx + x]) - row_mean;
    }

    row_rss = get_root_square_sum(&X[y*nx], nx);

    for (int x=0; x<nx; x++) {
      X[y*nx + x] = X[y*nx + x] / row_rss;
    }
  }


  // Allocate GPU memory for X
  CHECK_CUDA_ERROR( cudaMalloc( (void**)&dev_X, ny * nx * sizeof(double) ) );
  CHECK_CUDA_ERROR( cudaMalloc( (void**)&dev_result, ny * ny * sizeof(float) ) );

  // Copy X to the GPU memory.
  CHECK_CUDA_ERROR( cudaMemcpy(dev_X,
			       X,
			       sizeof(double)*nx*ny,
			       cudaMemcpyHostToDevice) );

  dim3 szBlock(8, 8);
  dim3 szGrid((ny + szBlock.x - 1) / szBlock.x,
  	      (ny + szBlock.y - 1) / szBlock.y);

  matrix_multiply<<<szGrid, szBlock>>>(nx, ny, dev_X, dev_result);

  // Copy result back to the CPU memory
  CHECK_CUDA_ERROR( cudaMemcpy(result,
  			       dev_result,
  			       sizeof(float)*ny*ny,
  			       cudaMemcpyDeviceToHost) );

  free(X);
  cudaFree(dev_X);
  cudaFree(dev_result);
}
