#include "so.h"
#include <omp.h>
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <iostream>

void merge_blocks(data_t* b1, data_t* b2, int n1, int n2) {
  data_t* b3 = (data_t*)malloc(sizeof(data_t) * n1);
  memcpy(b3, b1, sizeof(data_t)*n1);

  int i=0, j=0;
  for (int n=0; n<n1; n++) {
    if (b2[j] < b3[i]) {
      b1[n] = b2[j++];
    } else {
      b1[n] = b3[i++];
    }
  }
  for (int n=0; n<n2; n++) {
    if (i >= n1) {
      b2[n] = b2[j++];
    } else if (j >= n2) {
      b2[n] = b3[i++];
    } else {
      if (b2[j] < b3[i]) {
	b2[n] = b2[j++];
      } else {
	b2[n] = b3[i++];
      }
    }
  }

  free(b3);
}

void psort(int n, data_t* data) {
  if (n < 32) {
    std::sort(data, data + n);
    return;
  }

  // padding: the number of left over elements when dividing into p blocks
  int thread_count, blocksize, padding;

  int max_threads = omp_get_max_threads();

  // Start with power of 2 number of threads
  thread_count = std::pow(2, std::floor(std::log2(max_threads)));
  blocksize = std::floor(n / thread_count);
  padding = n - blocksize * thread_count;

  #pragma omp parallel num_threads(thread_count)
  {
  int thread_num = omp_get_thread_num();
  int len = (thread_num < thread_count-1) ? blocksize : blocksize + padding;
  std::sort(&data[thread_num*blocksize], &data[thread_num*blocksize] + len);
  }

  // Assuming that the thread count is a power of 2
  for (int p=thread_count; p > 1; p /= 2) {
    blocksize = std::floor(n / p);
    padding = n - blocksize * p;
    for (int i=0; i<p; i+=2) {
      int blocksize2 = (i < p-2) ? blocksize : blocksize + padding;//n1 + (n - blocksize1*p) : blocksize1;
      merge_blocks(&data[i*blocksize], &data[(i+1)*blocksize], blocksize, blocksize2);
    }
  }
}
