#include "cp.h"
#include <math.h>
#include <stdio.h>
#include <iostream>

#define POINT y*nx + x

double get_mean(const float* row, int nx) {
  double sum = 0.0, mean;
  for (int x=0; x<nx; x++) {
    sum += (double) row[x];
  }
  mean = sum / (double) nx;
  return mean;
}

double get_root_square_sum(const float* row, int nx) {
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
  double row_mean, row_rss, sum, square_sum, dp;
  double X[nx*ny], v1[nx], v2[ny];

  for (int y=0; y<ny; y++) {
    sum = 0.0;
    square_sum = 0.0;

    for (int x=0; x<nx; x++) {
      double paska = data[POINT];
      sum += paska;
    }

    row_mean = sum / ((double) nx);

    for (int x=0; x<nx; x++) {
      X[POINT] = ((double) data[POINT]) - row_mean;
      square_sum += X[POINT] * X[POINT];
    }

    row_rss = sqrt(square_sum);

    for (int x=0; x<nx; x++) {
      X[POINT] = X[POINT] / row_rss;
    }

  }

  for (int y=0; y<ny; y++) {
    for (int x=y; x<ny; x++) {
      for (int i=0; i<nx; i++) {
	v1[i] = X[y*ny + i];
      }

      for (int i=0; i<nx; i++) {
	v2[i] = X[x*ny + i];
      }

      dp = dot_product(v1, v2, nx);
      result[POINT] = (float) dp;
    }
  }
}
