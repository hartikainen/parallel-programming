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

void correlate(int ny, int nx, const float* data, float* result) {
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

#pragma omp parallel for schedule(static, 1)
  for (int y=0; y<ny; y++) {
    //printf("\n\nrow incoming: \n\n");
    for (int x=y; x<ny; x++) {
      double4_t r = {0.0, 0.0, 0.0, 0.0};
      double4_t* row1 = &X[y*nnx];
      double4_t* row2 = &X[x*nnx];

      for (int i=0; i<nnx; i++) {
	//r += X[x*nnx + i] * X[y*nnx + i];
	r += row1[i] * row2[i];

	// for (int x2=0; x2<4; x2++) {
	//   printf("row1: %f, row2: %f\n", row1[i][x2], row2[i][x2]);
	// }
      }
      double t = 0.0;
      for (int i=0; i<4; i++) {
	t += r[i];
      }
      result[y*ny + x] = t;
    }
  }
}
