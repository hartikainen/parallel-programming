#include "cp.h"
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
  double X[nx*ny];

  #pragma omp parallel for
  for (int y=0; y<ny; y++) {
    double row_mean = get_mean(&data[y*nx], nx);

    for (int x=0; x<nx; x++) {
      X[y*nx + x] = ((double) data[y*nx + x]) - row_mean;
    }

    double row_rss = get_root_square_sum(&X[y*nx], nx);

    for (int x=0; x<nx; x++) {
      X[y*nx + x] = X[y*nx + x] / row_rss;
    }
  }

  #pragma omp parallel for schedule(static, 1)
  for (int y=0; y<ny; y++) {
    for (int x=y; x<ny; x++) {
      double r = 0.0;
      for (int i=0; i<nx; i++) {
	r += X[x*nx + i] * X[y*nx + i];
      }
      result[y*ny + x] = r;
    }
  }
}
