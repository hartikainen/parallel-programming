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
  root_square_sum = sqrt(square_sum);
  return root_square_sum;
}

double dot_product(double* v1, double* v2, int len) {
  double dp = 0.0;

  for (int i=0; i<len; i++) {
    dp += v1[i] * v2[i];
  }

  return dp;
}

void correlate(int ny, int nx, const float* data, float* result) {
  double row_mean, row_rss, dp;
  double X[nx*ny], v1[nx], v2[nx];

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

  for (int y=0; y<ny; y++) {
    for (int x=y; x<ny; x++) {

      for (int i=0; i<nx; i++) {
	v1[i] = X[x*nx + i];
	v2[i] = X[y*nx + i];
      }

      dp = dot_product(v1, v2, nx);
      result[y*ny + x] = (float) dp;
    }
  }
}
