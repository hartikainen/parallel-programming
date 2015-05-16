#include "cp.h"
#include <cuda_runtime.h>
#include <cmath>
#include <stdio.h>
#include <iostream>

// input in global memory, get reasonable chunk to shared memory.
// shared memory -> registers -> do multiplications

// what are the bottlenecks?
// take kernel -> comment everything -> get the performance for different parts

// make -j parallelizes compilation
// make SMS=30: only produce machine code only from specific machine

// shared memory: 48kb per SM
// global memory: 2gb
// 2 SM's

#define CHECK_CUDA_ERROR(call) do { \
  cudaError_t result_ = (call); \
  if (result_ != cudaSuccess) { \
  fprintf(stderr, #call " failed: %s\n", \
	  cudaGetErrorString(result_)); \
  exit(1); \
  } \
} while(0)

float get_mean(const float* row, int nx) {
  float sum = 0.0, mean;
  for (int x=0; x<nx; x++) {
    sum += row[x];
  }
  mean = sum / (float) nx;
  return mean;
}

float get_root_square_sum(float* row, int nx) {
  float square_sum = 0, root_square_sum;
  for (int x=0; x<nx; x++) {
    square_sum += pow(row[x], 2.0);
  }
  root_square_sum = std::sqrt(square_sum);
  return root_square_sum;
}

__global__ void matrix_multiply(int nx, int ny, float* X, float* result, const int m, const int k) {
  if (x < y || x >= ny || y >= ny) return;

  __shared__ float A1[m][m][k];
  __shared__ float A2[m][m][k];
  float temp_X[k*k] = {0};
  float B1[k];
  float B2[k];

  float r = 0;

  int tx = threadIdx.x; // thread x in block
  int ty = threadIdx.y; // thread y in block
  int bx = blockIdx.x;
  int by = blockIdx.y;
  int x = tx + bx * blockDim.x; // global x
  int y = ty + by * blockDim.y; // global y

  for (int mx=0; mx<nx; mx+=m) {
    for (int i=0; i<k; i++) {
      // y =  k * m * blockIdx.x + threadIdx.y * k + i;
      // y =  k * m * blockIdx.y + threadIdx.y * k + i;
      A1[tx][ty][i] = X[k * m * bx + ty * k + i];
      A2[tx][ty][i] = X[k * m * by + ty * k + i];
    }

    __syncthreads();

    for (int i=0; i<k; i++) {
      for (int j=0; j<k; j++) {
	tmp[j*k + i] = s;
      }
    }
  }
  // for (int i=p; i<nx; i++) {
  //   r += X[x*nx + i] * X[y*nx + i];
  // }

  // result[y*ny + x] = r;
}

void correlate(int ny, int nx, const float* data, float* result) {
  float row_mean, row_rss;
  float* X = (float*) malloc(sizeof(float) * nx * ny);
  float* dev_X;
  float* dev_result;
  const int k = 4;
  const int m = 3;

  for (int y=0; y<ny; y++) {
    row_mean = get_mean(&data[y*nx], nx);

    for (int x=0; x<nx; x++) {
      X[y*nx + x] = (data[y*nx + x]) - row_mean;
    }

    row_rss = get_root_square_sum(&X[y*nx], nx);

    for (int x=0; x<nx; x++) {
      X[y*nx + x] = X[y*nx + x] / row_rss;
    }
  }


  // Allocate GPU memory for X
  CHECK_CUDA_ERROR( cudaMalloc( (void**)&dev_X, ny * nx * sizeof(float) ) );
  CHECK_CUDA_ERROR( cudaMalloc( (void**)&dev_result, ny * ny * sizeof(float) ) );

  // Copy X to the GPU memory.
  CHECK_CUDA_ERROR( cudaMemcpy(dev_X,
			       X,
			       sizeof(float)*nx*ny,
			       cudaMemcpyHostToDevice) );

  dim3 szBlock(m, m);
  dim3 szGrid((ny + szBlock.x - 1) / szBlock.x,
  	      (ny + szBlock.y - 1) / szBlock.y);

  matrix_multiply<<<szGrid, szBlock>>>(nx, ny, dev_X, dev_result, m, k);

  CHECK_CUDA_ERROR(cudaGetLastError());

  // Copy result back to the CPU memory
  CHECK_CUDA_ERROR( cudaMemcpy(result,
  			       dev_result,
  			       sizeof(float)*ny*ny,
  			       cudaMemcpyDeviceToHost) );

  free(X);
  cudaFree(dev_X);
  cudaFree(dev_result);
}
