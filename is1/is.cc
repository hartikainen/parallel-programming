#include "is.h"
#include <cmath>
#include <stdio.h>
#include <iostream>
#include "../common/vector.h"
#include <omp.h>

#define d4tod(d4) d4[0]+d4[1]+d4[2]

Result segment(int ny, int nx, const float* data) {
  int size = nx * ny; // size of the whole image
  double4_t* cdata = double4_alloc(size); // vectorized data
  double4_t vPc = {0.0};
  double4_t zero4_t = {0.0, 0.0, 0.0, 0.0};
  double maxmaxhXY = 0.0;

  int snx = nx+1, sny = ny+1;
  double4_t* S00 = double4_alloc(snx*sny);

  int ry0=0, rx0=0, ry1=1, rx1=1; // return coordinates
  double4_t rXc = {0.0}, rYc = {0.0}; // return colors
  Result result { ny/3, nx/3, 2*ny/3, 2*nx/3, {0.0f, 0.0f, 1.0f}, {1.0f, 0.0f, 0.0f} };

  #pragma omp for schedule(static, 1)
  for (int y=0; y<ny; y++) {
    for (int x=0; x<nx; x++) {
      int idx = y*nx + x;
      for (int c=0; c<3; c++) {
	cdata[idx][c] = data[c + 3*idx];
      }
    }
  }

  // The borders (x=0 and y=0) for the color sum matrix
  // are 0. Thus we can dynamically calculate the other sums.
  for (int sx=0; sx<snx; sx++) {
    S00[sx] = zero4_t;
  }
  for (int sy=0; sy<sny; sy++) {
    S00[snx*sy] = zero4_t;
  }
  for (int y1=0; y1<ny; y1++) {
    for (int x1=0; x1<nx; x1++) {
      S00[(y1+1)*snx + (x1+1)] = cdata[y1*nx + x1]
                               + S00[y1    *snx + (x1+1)]
                               + S00[(y1+1)*snx + x1    ]
                               - S00[y1    *snx + x1    ];
    }
  }

  vPc = S00[snx*sny-1];

#pragma omp parallel
{
  double4_t hXY4, vYc, vXc;
  double sizeX, sizeY, divX, divY, hXY, max_hXY = -1;
  int tx0 = 0, ty0 = 0, tx1 = 1, ty1 = 1;

  #pragma omp for schedule(static, 1)
  for (int h=1; h<=ny; h++) {
    for (int w=1; w<=nx; w++) {
      sizeX = (double)h * (double)w;
      sizeY = (double)size - sizeX;
      divX = 1.0/sizeX;
      divY = 1.0/sizeY;

      for (int y0=0; y0<=ny-h; y0++) {
	for (int x0=0; x0<=nx-w; x0++) {
	  int y1 = y0 + h;
	  int x1 = x0 + w;
	  vXc = S00[y1*snx + x1]
	      - S00[y1*snx + x0]
	      - S00[y0*snx + x1]
	      + S00[y0*snx + x0];
	  vYc = vPc - vXc;
	  hXY4 = (vXc * vXc * divX) + (vYc * vYc * divY);
	  hXY = d4tod(hXY4);

	  if (hXY > max_hXY) {
	    max_hXY = hXY;
	    tx0 = x0;
	    ty0 = y0;
	    tx1 = x1;
	    ty1 = y1;
	  }
	}
      }
    }
  }

  #pragma omp critical
  {
    if (max_hXY > maxmaxhXY) {
      maxmaxhXY = max_hXY;
      rx0 = tx0;
      ry0 = ty0;
      rx1 = tx1;
      ry1 = ty1;
    }
  }
}

#pragma omp parallel
{
  double4_t temp_rXc = {0.0}, temp_rYc = {0.0};
  #pragma omp for
  for (int y=0; y<ry0; y++) {
    for (int x=0; x<nx; x++) {
      temp_rYc += cdata[y*nx + x];
    }
  }
  #pragma omp for
  for (int y=ry0; y<ry1; y++) {
    for (int x=0; x<rx0; x++) {
      temp_rYc += cdata[y*nx + x];
    }
    for (int x=rx0; x<rx1; x++) {
      temp_rXc += cdata[y*nx + x];
    }
    for (int x=rx1; x<nx; x++) {
      temp_rYc += cdata[y*nx + x];
    }
  }
  #pragma omp for
  for (int y=ry1; y<ny; y++) {
    for (int x=0; x<nx; x++) {
      temp_rYc += cdata[y*nx + x];
    }
  }

  #pragma omp critical
  {
    rXc += temp_rXc;
    rYc += temp_rYc;
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
