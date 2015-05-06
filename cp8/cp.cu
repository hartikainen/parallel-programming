#include "cp.h"
#include <cuda_runtime.h>
#include <cmath>
#include <stdio.h>
#include <iostream>

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

void correlate(int ny, int nx, const float* data, float* result) {
  double row_mean, row_rss, X[nx*ny];

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

  double r, s;
  for (int y=0; y<ny; y++) {
    for (int x=y; x<ny; x++) {
      r=0;
      for (int i=0; i<nx; i++) {
	s = X[x*nx + i] * X[y*nx + i];
	r += s;
      }
      result[y*ny + x] = r;
    }
  }
}
