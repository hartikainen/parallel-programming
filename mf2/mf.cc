#include <algorithm>
#include "mf.h"

void mf(int ny, int nx, int hy, int hx, const float* in, float* out) {


#pragma omp parallel for
  // Loop through the pixels in the image
  for (py=0; py<ny; py++) {
    for (px=0; px<nx; px++) {

      /* px, py: loop indices for the actual points in the image
       * n: the index for assigning window values
       * i, j: loop indices for the window
       * wxmin, wxmax, wymin, wymax: the boundaries for the rectangle window
       * wsize: the size of the window array
       */
      int px, py, n, i, j, wxmin, wxmax, wymin, wymax, wsize;

      // window: The array to store the window values and count the median from
      float* window;
      // m1 and m2 used to store the intermediate median values in case that
      // the length of the median array is even
      float median, m1, m2;

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

      // Find the median. If the window's length is even, median
      // is the average of the two 'middle' values (m1, m2)
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

      out[py*nx + px] = median;
      free(window);
    }
  }
}
