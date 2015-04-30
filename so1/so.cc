#include "so.h"
#include <omp.h>
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <iostream>

// #pragma omp parallel num_threads(p)
// omp_get_max_threads()
// omp_get_thread_num()

void print_data(data_t* data, int n) {
  printf("\n data: \n");
  for (int i=0; i<n; i++) {
    printf("%.2f, ", (float)data[i]);
  }
  printf("\n");
}

void merge_blocks(data_t* b1, data_t* b2, int n1, int n2) {
  //  printf("\nmerging arrays: \n");
  //  for (int i=0; i<n1; i++) std::cout << b1[i] << " "; std::cout << "\n";
  // for (int i=0; i<n2; i++) std::cout << b2[i] << " "; std::cout << "\n";
  // printf("n1: %d, n2: %d\n", n1, n2);
  data_t* b3 = (data_t*)malloc(sizeof(data_t) * n1);
  memcpy(b3, b1, sizeof(data_t)*n1);
  // for (int n=0; n<n1; n++) {
  //   b3[n] = b1[n];
  // }

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
      // printf("\n\ni >= n1, overflow: n=%d, n1=%d, n2=%d, i=%d, j=%d, val=", n, n1, n2, i, j);
      // std::cout << b2[j] << "\n\n";
      b2[n] = b2[j++];
    } else if (j >= n2) {
      // printf("\n\nj >= n2, overflow: n=%d, n1=%d, n2=%d, i=%d, j=%d, val=", n, n1, n2, i, j);
      // std::cout << b2[j] << "\n\n";
      b2[n] = b3[i++];
    } else {
      if (b2[j] < b3[i]) {
	b2[n] = b2[j++];
      } else {
	b2[n] = b3[i++];
      }
    }
  }
  // printf("merge result:\n");
  // for (int iii=0; iii<n1+n2; iii++) std::cout << b1[iii] << " "; std::cout << "\n\n\n";
  free(b3);
}

void psort(int n, data_t* data) {
  data_t* data_copy = (data_t*)malloc(sizeof(data_t) * n);
  memcpy(data_copy, data, sizeof(data_t)*n);

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

  // printf("\nSORTEDBLOCKS:\n");
  // for (int j=0; j< thread_count; j++) {
  //   printf("[");
  //   int len = (j < thread_count-1) ? blocksize : blocksize + padding;
  //   for (int i=0; i<len; i++) std::cout << data[blocksize * j + i] << " "; std::cout << "";
  //   printf("], \n");
  // }


  // Assuming that the thread count is a power of 2
  for (int p=thread_count; p > 1; p /= 2) {
    blocksize = std::floor(n / p);
    padding = n - blocksize * p;
    // printf("P=%d, blocksize=%d, padding=%d\n", p, blocksize, padding);
    for (int i=0; i<p; i+=2) {
      int blocksize2 = (i < p-2) ? blocksize : blocksize + padding;//n1 + (n - blocksize1*p) : blocksize1;
      // printf("\ti=%d, blocksize2=%d\n", i, blocksize2);
      // printf("merging:\n");
      merge_blocks(&data[i*blocksize], &data[(i+1)*blocksize], blocksize, blocksize2);
    }
  }

  // std::sort(data_copy, data_copy+n);
  // for (int i=0; i<n; i++) {
  //   int d1 = data[i];
  //   int d2 = data_copy[i];
  //   std::cout << "got: " << d1 << ", should be: " << d2 << "\n";
  // }
}
