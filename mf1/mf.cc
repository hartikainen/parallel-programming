#include <algorithm>
#include "mf.h"

void mf(int ny, int nx, int hy, int hx, const float* in, float* out) {
  int px, py, n, i, j, wxmin, wxmax, wymin, wymax, wsize;
  float* window;
  float median, m1, m2;
  printf("nx*ny: %d\n", nx*ny);
  printf("hx:%d, hy:%d\n", hx, hy);
  // Loop through the pixels in the image
  for (py=0; py<ny; py++) {
    for (px=0; px<nx; px++) {
	// limit the window, s.t. it's smaller near the boundaries
      wxmin = (px - hx) > 0        ? (px - hx) : 0;
      wxmax = (px + hx) < (nx - 1) ? (px + hx+1) : nx;
      wymin = (py - hy) > 0        ? (py - hy) : 0;
      wymax = (py + hy) < (ny - 1) ? (py + hy+1) : ny;

      wsize = (wxmax - wxmin) * (wymax - wymin);
      window = (float*) malloc(sizeof(float) * wsize);

      // Store the values in window
      n = 0;
      for (i=wxmin; i<wxmax; i++) {
	for (j=wymin; j<wymax; j++) {
	  window[n++] = in[j*nx + i];
	}
      }

      // Find the median, if window's length is even, median
      // is the average of the two 'middle' values
      if (wsize % 2 == 0) {
	  std::nth_element(window, window + wsize/2, window + wsize);
	  m2 = window[wsize/2];
	  std::nth_element(window, (window + wsize/2)-1, window + wsize);
	  m1 = window[wsize/2-1];
	  median = (m1 + m2) / 2.0;
      } else {
	  std::nth_element(window, window + wsize/2, window + wsize);
	  median = window[wsize/2];
      }

      //      printf("px=%d, py=%d, wx=(%d,%d), wy=(%d,%d), median=%f\n", px, py, wxmin, wxmax, wymin, wymax, median);

      out[py*nx + px] = median;
      free(window);
    }
  }
}
