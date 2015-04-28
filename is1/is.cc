#include "is.h"
#include <cmath>
#include <stdio.h>
#include <iostream>
#include "../common/vector.h"

// For each pixel (x,y) and colour component c, we define the error error(y,x,c) as follows:
//   Let colour(y,x,c) = data[c + 3 * x + 3 * nx * y].
//   If (x,y) is located outside the rectangle: error(y,x,c) = outer[c] - colour(y,x,c).
//   If (x,y) is located inside the rectangle: error(y,x,c) = inner[c] - colour(y,x,c).
//   The total cost of the segmentation is the sum of squared errors, that is, the sum of error(y,x,c) * error(y,x,c) over all 0 <= c < 3 and 0 <= x < nx and 0 <= y < ny.


// TASK IS1 Your task is to find a segmentation that minimises the total cost.

// Implement an algorithm for image segmentation. You can use the naive algorithm that tries out all possible locations 0 <= y0 < y1 <= ny and 0 <= x0 < x1 <= nx for the rectangle and finds the best one.

// However, your implementation has to be reasonably efficient. In particular, make sure that for each possible location of the rectangle you only need to perform O(1) operations to evaluate how good this position is. That is, the total complexity for an n × n figure has to be O(n4). See the hints for suggestions.

// You are supposed to use the following technique to parallelise your code:

//   Multiple threads with OpenMP: Different threads try out different candidate locations.
//   Vector instructions: Store the colour image so that each pixel is a 4-element vector, with 3 colour components and a zero for padding. This way you can do all calculations for all colour components in parallel.

// Examples of a sufficient performance with classroom computers:
//   nx = ny = 400: roughly 10 seconds

#define d4tod(d4) d4[0]+d4[1]+d4[2]+d4[3]

void print_data(int ny, int nx, const float* data, double4_t* cdata) {
  printf("\n\nINITIAL DATA:\n");
  printf("nx: %d, ny: %d\n", nx, ny);
  for (int y=0; y<ny; y++) {
    for (int x=0; x<nx; x++) {
      printf("{");
      for (int c=0; c<3; c++) {
	printf("%f", data[c + 3 * x + 3 * nx * y]);
	if (c<2) printf(",");
      }
      printf("}");
    }
    printf("\n");
  }
  printf("\n");
  printf("\nVECTORIZED DATA\n");
  for (int y=0; y<ny; y++) {
    for (int x=0; x<nx; x++) {
      printf("{");
      for (int c=0; c<4; c++) {
	printf("%f", cdata[y*nx + x][c]);
	if (c<3) printf(",");
      }
      printf("}");
    }
    printf("\n");
  }
  printf("END OF DATA PRINT:\n");
}

void printS00(int ny, int nx, double* S00) {
  printf("\n\nS00:\n");
  printf("nx: %d, ny: %d\n", nx, ny);
  for (int y=0; y<ny; y++) {
    for (int x=0; x<nx; x++) {
      printf("%f ", S00[y*nx + x]);
    }
    printf("\n");
  }
  printf("\n");
  printf("END OF S00 PRINT:\n");
}

Result segment(int ny, int nx, const float* data) {
  int miny0, miny1, minx0, minx1;
  int size = nx * ny;
  double4_t cdata[ny*nx];
  double4_t vPc4 = {0};
  double vPc = 0;
  double vXc, vYc;
  double S00[ny*nx];
  double hXY = 0.0, temp_hXY;

  for (int y=0; y<ny; y++) {
    for (int x=0; x<nx; x++) {
      int idx = nx * y + x;
      for (int c=0; c<3; c++) {
	cdata[idx][c] = data[c + 3*idx];
      }
      vPc4 += cdata[idx];
    }
  }

  vPc = d4tod(vPc4);

  // TODO: remove
  print_data(ny, nx, data, cdata);
  printf("vPc=%f", vPc);


  // Initialize the borders (x=0 and y=0) for the
  // color sum matrix. Then we can dynamically calculate
  // the other sums.
  S00[0] = d4tod(cdata[0]);
  for (int x1=0; x1<nx; x1++) {
    S00[x1] = d4tod(cdata[x1]) + S00[x1-1];
  }
  for (int y1=0; y1<ny; y1++) {
    S00[y1*nx] = d4tod(cdata[y1*nx]) + S00[(y1-1)*nx];
  }

  for (int y1=1; y1<ny; y1++) {
    for (int x1=1; x1<nx; x1++) {
      S00[y1*nx + x1] = d4tod(cdata[y1*nx + x1]) + S00[(y1-1)*nx + (x1)] + S00[y1*nx + (x1-1)] - S00[(y1-1)*nx + (x1-1)];
    }
  }

  // TODO: remove
  printS00(ny, nx, S00);

  for (int h=1; h<=ny; h++) {
    for (int w=1; w<=nx; w++) {
      printf("h=%d, w=%d: ", h, w);
      int sizeX = h * w;
      int sizeY = size - sizeX;

      for (int y0=0; y0<=ny-h; y0++) {
	for (int x0=0; x0<=nx-w; x0++) {
	  int y1 = y0 + h;
	  int x1 = x0 + w;
	  vXc = S00[y1*nx + x1]
	      - S00[y1*nx + x0]
	      - S00[y0*nx + x1]
	      + S00[y0*nx + x0];
	  vYc = vPc - vXc;
	  temp_hXY = ((vXc*vXc) / sizeX) + ((vYc*vYc) / sizeY);
	  // printf(" temp_hXY: %f", temp_hXY);
	  printf("y0: %d, x0: %d, y1: %d, x1: %d, hXy: %f, vXc: %f, vYc: %f --  ", y0, x0, y1, x1, temp_hXY, vXc, vYc);

	  if (temp_hXY > hXY) {
	    hXY = temp_hXY;
	    miny0 = y0; miny1 = y1; minx0 = x0; minx1 = x1;
	    // printf("\ninside if\n");
	    // printf("x0: %d, y0: %d, w: %d, h: %d, hXY: %f\n", x0, y0, w, h, hXY);
	  }
	}
      }
    }
    printf("\n");
  }

  printf("RESULT -> minx0: %d, miny0: %d, minx1: %d, miny1: %d, minhXY: %f\n", minx0, miny0, minx1, miny1, hXY);


  Result result { ny/3, nx/3, 2*ny/3, 2*nx/3, {0.0f, 0.0f, 1.0f}, {1.0f, 0.0f, 0.0f} };
  return result;
}
