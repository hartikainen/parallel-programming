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

#define M 2
#define K 2
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
  int x = tx + bx * blockDim.x; // global x
  int y = ty + by * blockDim.y; // global y

  int x_base = m*k*bx + ty*k; // y1
  int y_base = m*k*by + ty*k; // y2


  __syncthreads();
  if (tx == 0 && ty == 0) printf("block at: y=%d, x=%d\n", by, bx);
  __syncthreads();
  printf("thread at point: y=%d, x=%d\n", y, x);
  __syncthreads();
  printf("txty: ty=%d, tx=%d\n", ty, tx);
  __syncthreads();

  // if (x < y || x >= ny || y >= ny) return;

  for (int mx=0; mx<nx; mx+=m) {
    for (int i=0; i<k; i++) {
      // Are these correct??
      A1[ty][tx][i] = X[(x_base + i)*nx + (mx + tx)];
      A2[ty][tx][i] = X[(y_base + i)*nx + (mx + tx)];
    }

    if (tx == 0 && ty == 0) printf("printing cache, mx=%d, tx=%d, ty=%d", mx, tx, ty);
    for (int j=0; j<m; j++) {
      for (int kk=0; kk<k; kk++) {
	printf("\n");
	for (int i=0; i<m; i++) {
	  if (tx == 0 && ty == 0) printf("%f\t", A1[j][i][kk]);
	}
      }
      if (tx == 0 && ty == 0) printf("\n");
    }

    __syncthreads();

    for (int mm=0; mm<m; mm++) {
      if (tx == 0 && ty == 0) printf("mm=%d\n", mm);
      for (int i=0; i<k; i++) {
	for (int j=0; j<k; j++) {
	  tmp[j][i] += A1[tx][mm][i] * A2[ty][mm][j];
	  if (tx == 0 && ty == 0) printf("%f\t", tmp[j][i]);
	}
	if (tx == 0 && ty == 0) printf("\n");
      }
    }
  }

  for (int i=0; i<k; i++) {
    for (int j=0; j<k; j++) {
      result[(y_base + j) * ny + (x_base + i)] = tmp[j][i];
    }
  }
}

void correlate(int ny, int nx, const float* data, float* result) {
  float row_mean, row_rss;
  float* X = (float*) malloc(sizeof(float) * nx * ny);
  float* dev_X;
  float* dev_result;
  int m = M;
  int k = K;

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

  printf("\nnormalized data:\n");
  for (int j=0; j<ny; j++) {
    for (int i=0; i<nx; i++) {
      std::cout << X[j*nx + i] << " ";
    }
    printf("\n");
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
  dim3 szGrid((ny + (szBlock.x*k) - 1) / (szBlock.x * k),
  	      (ny + (szBlock.y*k) - 1) / (szBlock.y * k));

  matrix_multiply<<<szGrid, szBlock>>>(nx, ny, dev_X, dev_result);

  CHECK_CUDA_ERROR(cudaGetLastError());

  // Copy result back to the CPU memory
  CHECK_CUDA_ERROR( cudaMemcpy(result,
  			       dev_result,
  			       sizeof(float)*ny*ny,
  			       cudaMemcpyDeviceToHost) );


  printf("\nresult:\n");
  for (int j=0; j<ny; j++) {
    for (int i=0; i<ny; i++) {
      std::cout << result[j*nx + i] << " ";
    }
    printf("\n");
  }

  free(X);
  cudaFree(dev_X);
  cudaFree(dev_result);
}
