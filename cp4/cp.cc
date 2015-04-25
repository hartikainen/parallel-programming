#include "cp.h"
#include <cmath>
#include <stdio.h>
#include <iostream>
#include "../common/vector.h"

#define BLOCK_SIZE 2

double get_sum(double4_t* row, int n) {
  double4_t sum4_t = {0.0, 0.0, 0.0, 0.0};
  double sum = 0.0;

  for (int i=0; i<n; i++) {
    sum4_t += row[i];
  }

  for (int i=0; i<4; i++) {
    sum += sum4_t[i];
  }

  return sum;
}

double get_root_square_sum(double4_t* row, int n) {
  double row_ss = 0.0, row_rss;
  double4_t tmp4, row_ss4 = {0.0, 0.0, 0.0, 0.0};

  for (int i=0; i<n; i++) {
    tmp4 = row[i];
    row_ss4 += tmp4 * tmp4;
  }

  for (int i=0; i<4; i++) {
    row_ss += row_ss4[i];
  }

  row_rss = std::sqrt(row_ss);

  return row_rss;
}

double dot_product(double4_t* v1, double4_t* v2, int len) {
  double sum = 0.0;
  double4_t sum4 = {0.0, 0.0, 0.0, 0.0};

  for (int i=0; i<len; i++) {
    sum4 += v1[i] * v2[i];
  }

  for (int i=0; i<4; i++) {
    sum += sum4[i];
  }

  return sum;
}

void block_dot_product(double4_t* B1[], double4_t* B2[], int bs, int len, float* result) {
  double4_t tmp[bs*bs] = {{0.0}};
  for (int i=0; i<len; i++) {
    for (int j1=0; j1<bs; j1++) {
      for (int j2=0; j2<bs; j2++) {
	tmp[j1*bs + j2] += B1[j1][i] * B2[j2][i];
      }
    }
  }

  for (int i=0; i<4; i++) {
    for (int j=0; j<bs; j++) {
      //B3[j] += tmp[j][i];
      result[(y+j)*ny + (x+i)] += tmp[j][i];//B3[j*bs + i];
    }
  }
}

void correlate(int ny, int nx, const float* data, float* result) {
  int bs = BLOCK_SIZE;
  int nnx = std::ceil(nx/4.0);
  int xpad = nnx * 4.0 - nx;
  double4_t* X = double4_alloc(nnx*ny);

  // Copy the data into double4_t array X
#pragma omp parallel for
  for (int y=0; y<ny; y++) {
    for (int x1=0; x1<nnx-1; x1++) {
      for (int x2=0; x2<4; x2++) {
	X[y*nnx + x1][x2] = data[y*nx + x1*4 + x2];
      }
    }

    for (int x=0; x<4-xpad; x++) {
      X[y*nnx + nnx-1][x] = data[y*nx + (nnx - 1) * 4 + x];
    }
    // Initialize the remaining 'right' padding to 0.0
    for (int x=4-xpad; x<4; x++) {
      X[y*nnx + nnx-1][x] = 0.0;
    }
  }

#pragma omp parallel for
  for (int y=0; y<ny; y++) {
    double row_sum, row_mean, row_rss;

    row_sum = get_sum(&X[y*nnx], nnx);
    row_mean = row_sum / (double)nx;

    for (int x1=0; x1<nnx-1; x1++) {
      X[y*nnx + x1] = X[y*nnx + x1] - row_mean;
    }


    for (int x=0; x<4-xpad; x++) {
      X[y*nnx + nnx-1][x] = data[y*nx + (nnx - 1) * 4 + x] - row_mean;
    }

    row_rss = get_root_square_sum(&X[y*nnx], nnx);

    for (int x=0; x<nnx; x++) {
      X[y*nnx + x] /= row_rss;
    }
  }

  int boundary = std::floor(ny/bs)*bs;
#pragma omp parallel for schedule(static, 1)
  for (int y=0; y<boundary; y+=bs) {
    for (int x=y; x<boundary; x+=bs) {
      double4_t* B1[bs];
      double4_t* B2[bs];

      for (int i=0; i<bs; i++) {
	B1[i] = &X[(y+i)*nnx];
	B2[i] = &X[(x+i)*nnx];
      }

      //  double B3[bs*bs] = {0.0};
      block_dot_product(B1, B2, bs, nnx, &result[y*ny + x]);

      for (int j=0; j<bs; j++) {
	for (int i=0; i<bs; i++) {

	}
      }
    }
  }

  for (int j=0; j<ny; j++) {
    double4_t* v1 = &X[j*nnx];
    double4_t* v2;
    for (int i=boundary; i<ny; i++) {
      v2 = &X[i*nnx];
      result[j*ny + i] = dot_product(v1, v2, nnx);
    }
  }
}
