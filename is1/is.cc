#include "is.h"
#include <cmath>
#include <stdio.h>
#include <iostream>
#include "../common/vector.h"

#define d4tod(d4) d4[0]+d4[1]+d4[2]+d4[3]

Result segment(int ny, int nx, const float* data) {
  int size = nx * ny; // size of the whole image
  double4_t cdata[size] = {0.0}; // vectorized data
  double4_t vPc4 = {0.0};

  double vPc = 0.0, max_hXY = 0.0;
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
      double sizeX = (double)h * (double)w;
      double sizeY = (double)size - sizeX;
      double divX = 1.0/sizeX;
      double divY = 1.0/sizeY;

      for (int y0=0; y0<=ny-h; y0++) {
	for (int x0=0; x0<=nx-w; x0++) {
	  int y1 = y0 + h;
	  int x1 = x0 + w;
	  double hXY, vYc, vXc;
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
  rXc /= (double)sizeX;
  rYc /= (double)(size - sizeX);

  result.y0 = ry0;
  result.x0 = rx0;
  result.x1 = rx1;
  result.y1 = ry1;
  for (int c=0; c<3; c++) {
    result.inner[c] = (float)rXc[c];
    result.outer[c] = (float)rYc[c];
  }

  return result;
}
