#include "cp.h"
#include <cmath>
#include <stdio.h>
#include <iostream>
#include "../common/vector.h"

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

double4_t get_root_square_sum(double4_t* row, int n) {
  double row_ss = 0.0, row_rss;
  double4_t tmp4, row_ss4 = {0.0};

  for (int i=0; i<n; i++) {
    tmp4 = row[i];
    row_ss4 += tmp4 * tmp4;
  }

  for (int i=0; i<4; i++) {
    row_ss += row_ss4[i];
  }

  row_rss = std::sqrt(row_ss);

  double4_t row_rss4 = {row_rss, row_rss, row_rss, row_rss};
  return row_rss4;
}

void correlate(int ny, int nx, const float* data, float* result) {
  int nnx = std::ceil(nx/4.0);
  int xpad = nnx * 4.0 - nx;
  double4_t* X = double4_alloc(nnx*ny);

  // Copy the data into double4_t array X
  for (int y=0; y<ny; y++) {
    for (int x1=0; x1<nnx-1; x1++) {
      for (int x2=0; x2<4; x2++) {
	X[y*nnx + x1][x2] = data[y*nx + x1*4 + x2];
	printf("times here: %d", x2);
      }
    }

    for (int x=0; x<xpad; x++) {
      X[y*nnx + nnx-1][x] = data[y*nx + (nnx - 1) * 4 + x];
    }
    // Initialize the remaining 'right' padding to 0.0
    for (int x=4-xpad; x<4; x++) {
      X[y*nnx + nnx-1][x] = 0.0;
    }
  }

  // std::cout << "\ndata: \n";
  // for (int y=0; y<ny; y++) {
  // 	for (int x=0; x<nx; x++) {
  // 	    printf("%f ", data[y*nx + x]);
  // 	}
  // 	printf("\n");
  // }

  // std::cout << "\nX: \n";
  // for (int y=0; y<ny; y++) {
  // 	for (int x1=0; x1<nnx; x1++) {
  // 	    for (int x2=0; x2<4; x2++) {
  // 		printf("%f ", X[y*nnx + x1][x2]);
  // 	    }
  // 	}
  // 	printf("\n");
  // }

#pragma omp parallel for
  for (int y=0; y<ny; y++) {
    double row_sum, row_mean;
    double4_t row_rss;

    row_sum = get_sum(&X[y*nnx], nnx);
    row_mean = row_sum / (double)nx;

    double4_t row_mean4 = {row_mean, row_mean, row_mean, row_mean};;

    for (int x1=0; x1<nnx-1; x1++) {
      for (int x2=0; x2<4; x2++) {
	X[y*nnx + x1] = X[y*nnx + x1] - row_mean4;
      }
    }

    for (int x=0; x<xpad; x++) {
      X[y*nnx + nnx-1][x] = data[y*nx + (nnx - 1) * 4 + x] - row_mean;
    }

    row_rss = get_root_square_sum(&X[y*nnx], nnx);

    for (int x=0; x<nnx; x++) {
      X[y*nnx + x] /= row_rss;
    }
  }

  // std::cout << "\nzero mean X: \n";
  // for (int y=0; y<ny; y++) {
  //   for (int x1=0; x1<nnx; x1++) {
  //     for (int x2=0; x2<4; x2++) {
  // 	printf("%f ", X[y*nnx + x1][x2]);
  //     }
  //   }
  //   printf("\n");
  // }

// #pragma omp parallel for schedule(static, 1)
//   for (int y=0; y<ny; y++) {
//     for (int x=y; x<ny; x++) {
//       double r = 0.0;
//       for (int i=0; i<nx; i++) {
// 	r += X[x*nx + i] * X[y*nx + i];
//       }
//       result[y*ny + x] = r;
//     }
//   }
}
