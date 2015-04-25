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

void dot_product2(double4_t* v1, double4_t* v2, double4_t* v3, double4_t* v4, int len, double* sum) {
  double4_t sum1 = {0.0, 0.0, 0.0, 0.0}, sum2 = {0.0, 0.0, 0.0, 0.0}, sum3 = {0.0, 0.0, 0.0, 0.0}, sum4 = {0.0, 0.0, 0.0, 0.0};

  // printf("v1: %f, %f, %f, %f\n", v1[1][0], v1[1][1], v1[1][2], v1[1][3]);
  // printf("v2: %f, %f, %f, %f\n", v2[1][0], v2[1][1], v2[1][2], v2[1][3]);
  for (int i=0; i<len; i++) {
    sum1 += v1[i] * v3[i];
    sum2 += v1[i] * v4[i];
    sum3 += v2[i] * v3[i];
    sum4 += v2[i] * v4[i];
  }

  for (int i=0; i<4; i++) {
    sum[0] += sum1[i];
    sum[1] += sum2[i];
    sum[2] += sum3[i];
    sum[3] += sum4[i];
  }

  // printf("jiiri:\n");
  // for (int i=0; i<4; i++) {printf("%f ", sum[i]);}
  // printf("\n");
}

void correlate(int ny, int nx, const float* data, float* result) {
  int bs = BLOCK_SIZE;
  int nnx = std::ceil(nx/4.0);
  int xpad = nnx * 4.0 - nx;
  double4_t* X = double4_alloc(nnx*ny);

  // std::cout << "\ndata: \n";
  // for (int y=0; y<ny; y++) {
  //   for (int x=0; x<nx; x++) {
  //     printf("%f ", data[y*nx + x]);
  //   }
  // }
  //   printf("\n");

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

  int asdf = std::floor(ny/bs)*bs;
#pragma omp parallel for schedule(static, 1)
  for (int y=0; y<asdf; y+=bs) {
    for (int x=y; x<asdf; x+=bs) {

      double4_t* v1 = &X[y*nnx];
      double4_t* v2 = &X[(y+1)*nnx];
      double4_t* v3 = &X[x*nnx];
      double4_t* v4 = &X[(x+1)*nnx];

      double jiiri[4] = {0.0, 0.0, 0.0, 0.0};
      dot_product2(v1, v2, v3, v4, nnx, jiiri);

      result[y*ny + x] = jiiri[0];
      result[y*ny + x+1] = jiiri[1];
      result[(y+1)*ny + x] = jiiri[2];
      result[(y+1)*ny + x+1] = jiiri[3];


      // for (int i=asdf; i<ny; i++) {
      // 	result[y*ny + i] = dot_product(&X[y*nnx], &X[i*nnx], nnx);
      // }

      // double4_t* v5;
      // double4_t* v6;
      // double4_t* v7;
      // double4_t* v8;
      // double4_t* v9;
      // double4_t r = {0.0, 0.0, 0.0, 0.0};
      //      int sz;
      //double t = 0.0;
    }
  }

  for (int j=0; j<ny; j++) {
    double4_t* v1 = &X[j*nnx];
    double4_t* v2;
    for (int i=std::floor(ny/bs)*bs; i<ny; i++) {
      v2 = &X[i*nnx];
      result[j*ny + i] = dot_product(v1, v2, nnx);
    }
  }

  // printf("result:\n");
  // for (int j=0; j<ny; j++) {
  //   for (int i=0; i<j; i++) {printf("........ ");}
  //   for (int i=j; i<ny; i++) {
  //     printf("%f ", result[j*ny + i]);
  //   }
  //   printf("\n");
  // }

}
