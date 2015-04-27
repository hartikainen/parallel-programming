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

Result segment(int ny, int nx, const float* data) {
  double4_t cdata[ny*nx];

  for (int y=0; y<ny; y++) {
    for (int x=0; x<nx; x++) {
      for (int c=0; c<3; c++) {
	cdata[y*nx + x][c] = data[c + 3 * x + 3 * nx * y];
      }
    }
  }

  print_data(ny, nx, data, cdata);

  Result result { ny/3, nx/3, 2*ny/3, 2*nx/3, {0.0f, 0.0f, 1.0f}, {1.0f, 0.0f, 0.0f} };
  return result;
}
