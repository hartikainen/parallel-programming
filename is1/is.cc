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

Result segment(int ny, int nx, const float* data) {
  int size = nx * ny; // size of the whole image
  double4_t cdata[size] = {0.0}; // vectorized data
  double4_t vPc4 = {0.0};

  double vXc, vYc, hXY, vPc = 0.0, max_hXY = 0.0;
  double S00[(ny+1)*(nx+1)] = {0.0};

  int ry0=0, rx0=0, ry1=1, rx1=1; // return coordinates
  double4_t rXc = {0.0}, rYc = {0.0}; // return colors
  Result result { ny/3, nx/3, 2*ny/3, 2*nx/3, {0.0f, 0.0f, 1.0f}, {1.0f, 0.0f, 0.0f} };

  for (int y=0; y<ny; y++) {
    for (int x=0; x<nx; x++) {
      int idx = y*nx + x;
      for (int c=0; c<3; c++) {
	cdata[idx][c] = data[c + 3*idx];
      }
      vPc4 += cdata[idx];
    }
  }

  vPc = d4tod(vPc4);

  // The borders (x=0 and y=0) for the color sum matrix
  // are 0. Thus we can dynamically calculate the other sums.
  for (int y1=0; y1<ny; y1++) {
    for (int x1=0; x1<nx; x1++) {
      S00[(y1+1)*(nx+1) + (x1+1)] = d4tod(cdata[y1*nx + x1])
                                  + S00[y1*(nx+1) + (x1+1)]
 	                          + S00[(y1+1)*(nx+1) + x1]
                                  - S00[y1*(nx+1) + x1];
    }
  }

  //#pragma omp parallel for schedule(static, 1)
  for (int h=1; h<=ny; h++) {
    for (int w=1; w<=nx; w++) {
      int sizeX = h * w;
      int sizeY = size - sizeX;
      double divX = 1.0/sizeX;
      double divY = 1.0/sizeY;

      for (int y0=0; y0<=ny-h; y0++) {
	for (int x0=0; x0<=nx-w; x0++) {
	  int y1 = y0 + h;
	  int x1 = x0 + w;
	  vXc = S00[y1*(nx+1) + x1]
	      - S00[y1*(nx+1) + x0]
	      - S00[y0*(nx+1) + x1]
	      + S00[y0*(nx+1) + x0];
	  vYc = vPc - vXc;
	  double tmpX = vXc * vXc * divX;
	  double tmpY = vYc * vYc * divY;
	  hXY = tmpX + tmpY;

	  if (hXY > max_hXY) {
	    max_hXY = hXY;
	    rx0 = x0;
	    ry0 = y0;
	    rx1 = x1;
	    ry1 = y1;
	  }
	}
      }
    }
  }

  for (int y=0; y<ry0; y++) {
    for (int x=0; x<nx; x++) {
      rYc += cdata[y*nx + x];
    }
  }
  for (int y=ry0; y<ry1; y++) {
    for (int x=0; x<rx0; x++) {
      rYc += cdata[y*nx + x];
    }
    for (int x=rx0; x<rx1; x++) {
      rXc += cdata[y*nx + x];
    }
    for (int x=rx1; x<nx; x++) {
      rYc += cdata[y*nx + x];
    }
  }

  for (int y=ry1; y<ny; y++) {
    for (int x=0; x<nx; x++) {
      rYc += cdata[y*nx + x];
    }
  }

  int sizeX = (ry1-ry0) * (rx1-rx0);
  rXc /= sizeX;
  rYc /= size - sizeX;

  result.y0 = ry0;
  result.x0 = rx0;
  result.x1 = rx1;
  result.y1 = ry1;
  for (int c=0; c<3; c++) {
    result.inner[c] = rXc[c];
    result.outer[c] = rYc[c];
  }

  return result;
}
