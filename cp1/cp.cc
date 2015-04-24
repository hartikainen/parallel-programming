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
  std::cout << "data: " << data[0] << " " << data[1] << " " << data[2] << " " << data[3] << "\n";
  std::cout << "nx: " << nx << " ny: " << ny << "\n";
  for (int y=0; y<ny; y++) {
    sum = 0.0;
    square_sum = 0.0;

    for (int x=0; x<nx; x++) {
      sum += (double) data[POINT];
    }

    std::cout << "\nrow sum=" << sum;
    row_mean = sum / ((double) nx);

    std::cout << "\nrow_mean=" << row_mean << "\n";

    for (int x=0; x<nx; x++) {
      X[POINT] = ((double) data[POINT]) - row_mean;
      square_sum += X[POINT] * X[POINT];
    }

    std::cout << "X after zero meaning:\n";
    std::cout << "X: " << X[0] << " " << X[1] << " " << X[2] << " " << X[3] << "\n";

    row_rss = sqrt(square_sum);

    std::cout << "row_rss=" << row_rss << "\n";

    for (int x=0; x<nx; x++) {
      X[POINT] = X[POINT] / row_rss;
    }

    std::cout << "X after unit variancing:\n";
    std::cout << "X: " << X[0] << " " << X[1] << " " << X[2] << " " << X[3] << "\n";

  }
  std::cout << "X: " << X[0] << " " << X[1] << " " << X[2] << " " << X[3] << "\n";

  printf("X:\n");
  for (int y=0; y<ny; y++) {
    //    for (int x=0; x<y; x++) printf("%.8f ", X[POINT]);
    for (int x=y; x<ny; x++) {
      printf("%.8f ", X[POINT]);

      for (int i=0; i<nx; i++) {
	v1[i] = X[y*ny + i];
      }

      for (int i=0; i<nx; i++) {
	v2[i] = X[x*ny + i];
      }

      dp = dot_product(v1, v2, nx);
      result[POINT] = (float) dp;
    }

    printf("\n");
  }

  printf("result:\n");
  for (int y=0; y<ny; y++) {
    for (int x=0; x<y; x++) printf(".......... ");
    for (int x=y; x<nx; x++) {
      printf("%.8f ", result[POINT]);
    }
    printf("\n");
  }
}
