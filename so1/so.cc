#include "so.h"
#include <omp.h>
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <iostream>

void merge_blocks(data_t* b1, data_t* b2, int n1, int n2) {
  data_t* b3 = (data_t*)malloc(sizeof(data_t) * n1);
  int max_threads = omp_get_max_threads();
  int part_size = n1 / max_threads;

  #pragma omp parallel for num_threads(max_threads)
  for (int t=0; t<max_threads; t++) {
    int nn = (t == max_threads) ? part_size : part_size + (n1 - part_size * max_threads);
    memcpy(&b3[t*part_size], &b1[t*part_size], sizeof(data_t)*nn);
  }


  int i=0, j=0, n=0;
  while (i < n1 && j < n2) {
    b1[n++] = (b2[j] < b3[i]) ? b2[j++] : b3[i++];
  }
  while (i < n1) {
    b1[n++] = b3[i++];
  }
  while (j < n2) {
    b1[n++] = b2[j++];
  }

  free(b3);
}

void psort(int n, data_t* data) {
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

    #pragma omp parallel for num_threads(p/2)
    for (int i=0; i<p; i+=2) {
      int blocksize2 = (i < p-2) ? blocksize : blocksize + padding;
      merge_blocks(&data[i*blocksize], &data[(i+1)*blocksize], blocksize, blocksize2);
    }
  }
}
