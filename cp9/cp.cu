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

#define M 25
#define K 5
__global__ void matrix_multiply(int nx, int ny, float* X, float* result) {
  const int m = M;
  const int k = K;
  __shared__ float A1[m][m][k];
  __shared__ float A2[m][m][k];
  float tmp[k][k] = {0};

  int tx = threadIdx.x; // thread x in block
  int ty = threadIdx.y;
  int bx = blockIdx.x;
  int by = blockIdx.y;

  int x_base = m*k*bx;
  int y_base = m*k*by;

  if (bx < by) return;

  for (int mx=0; mx<nx; mx+=m) {
    for (int i=0; i<k; i++) {
      A1[ty][tx][i] = X[(x_base + ty*k + i)*nx + (mx + tx)];
      A2[ty][tx][i] = X[(y_base + ty*k + i)*nx + (mx + tx)];
    }

    __syncthreads();

    for (int mm=0; mm<m; mm++) {
      for (int i=0; i<k; i++) {
	for (int j=0; j<k; j++) {
	  tmp[j][i] += A1[tx][mm][i] * A2[ty][mm][j];
	}
      }
    }
  }

  int y_base2 = y_base + ty*k;
  int x_base2 = x_base + tx*k;
  for (int i=0; i<k; i++) {
    int x_ind = x_base2 + i;
    for (int j=0; j<k; j++) {
      int y_ind = y_base2 + j;
      result[y_ind * ny + x_ind] = tmp[j][i];
    }
  }
}

void correlate(int ny, int nx, const float* data, float* result) {
  float row_mean, row_rss;
  float* dev_X;
  float* dev_result;
  int m = M;
  int k = K;

  int new_x = (nx/(m*k) + 1) * (m*k);
  int new_y = (ny/(m*k) + 1) * (m*k);

  dim3 szBlock(m, m);
  dim3 szGrid((new_y + (szBlock.x*k) - 1) / (szBlock.x * k),
  	      (new_y + (szBlock.y*k) - 1) / (szBlock.y * k));

  float* X = (float*) malloc(sizeof(float) * new_x * new_y);

  for (int y=0; y<ny; y++) {
    row_mean = get_mean(&data[y*nx], nx);

    for (int x=0; x<nx; x++) {
      X[y*new_x + x] = (data[y*nx + x]) - row_mean;
    }
    // pad the columns
    for (int x=nx; x<new_x; x++) {
      X[y*new_x + x] = 0.0;
    }

    row_rss = get_root_square_sum(&X[y*new_x], nx);

    for (int x=0; x<nx; x++) {
      X[y*new_x + x] = X[y*new_x + x] / row_rss;
    }
  }

  // pad the rows3
  for (int y=ny; y<new_y; y++) {
    for (int x=0; x<new_x; x++) {
      X[y*new_x + x] = 0.0;
    }
  }

  // Allocate GPU memory for X
  CHECK_CUDA_ERROR( cudaMalloc( (void**)&dev_X, new_y * new_x * sizeof(float) ) );
  CHECK_CUDA_ERROR( cudaMalloc( (void**)&dev_result, new_y * new_y * sizeof(float) ) );

  // Copy X to the GPU memory.
  CHECK_CUDA_ERROR( cudaMemcpy(dev_X,
			       X,
			       sizeof(float)*new_x*new_y,
			       cudaMemcpyHostToDevice) );

  matrix_multiply<<<szGrid, szBlock>>>(new_x, new_y, dev_X, dev_result);

  CHECK_CUDA_ERROR(cudaGetLastError());

  // Copy result back to the CPU memory
  for (int y=0; y<ny; y++) {
    CHECK_CUDA_ERROR( cudaMemcpy(&result[y*ny],
  				 &dev_result[y*new_y],
  				 sizeof(float)*ny,
  				 cudaMemcpyDeviceToHost) );
  }

  free(X);
  cudaFree(dev_X);
  cudaFree(dev_result);
}
