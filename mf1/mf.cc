#include <algorithm>
#include "mf.h"

void mf(int ny, int nx, int hy, int hx, const float* in, float* out) {
  int px, py, n, i, j, wxmin, wxmax, wymin, wymax, wsize;
  float* window;
  float median;

  // Loop through the pixels in the image
  for (py=0; py<ny; py++) {
    for (px=0; px<nx; px++) {
      wxmin = (px - hx) > 0        ? (px - hx) : 0;
      wxmax = (px + hx) < (nx - 1) ? (px + hx) : nx - 1;
      wymin = (py - hy) > 0        ? (py - hy) : 0;
      wymax = (py + hy) < (ny - 1) ? (py + hy) : ny - 1;
	  
      wsize = (wxmax - wxmin + 1) * (wymax - wymin + 1);
      window = (float*) malloc(wsize);

      n = 0;
      for (i=wxmin; i<wxmax+1; i++) {
	for (j=wymin; j<wymax+1; j++) {
	  window[n++] = in[j*nx + i];
	}
      }

      std::nth_element(window, window + wsize/2, window + wsize);
      median = window[wsize/2];

      out[py*nx + px] = median;
      free(window);
    }
  }
}
