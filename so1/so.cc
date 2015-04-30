#include "so.h"
#include <omp.h>
#include <algorithm>
#include <stdio.h>
#include <iostream>

// #pragma omp parallel num_threads(p)
// omp_get_max_threads()
// omp_get_thread_num()

void print_actual_data(data_t* data, int n) {
  printf("\n\n actual data: \n");
  for (int i=0; i<n; i++) {
    printf("%10.2f, ", (float)data[i]);
  }

  printf("\n");
}

void print_blocks(data_t** blocks, int bs, int p, int n) {
  printf("\n");
  printf("blocks: \n");
  for (int j=0; j<p-1; j++) {
    for (int i=0; i<bs; i++) {
      printf("%10.2f, ", (float)blocks[j][i]);
    }
    printf("\n\n");
  }
  for (int i=0; i<bs + (n - bs * p); i++) {
    printf("%10.2f, ", (float)blocks[p-1][i]);
  }
  printf("\n");
}

void merge_blocks(data_t* b1, data_t* b2, int n1, int n2) {
  printf("\nmerging arrays: \n");
  for (int i=0; i<n1; i++) std::cout << b1[i] << " "; std::cout << "\n";
  for (int i=0; i<n2; i++) std::cout << b2[i] << " "; std::cout << "\n\n\n";
  data_t b3[n1];
  for (int n=0; n<n1; n++) {
    b3[n] = b1[n];
  }

  int i=0, j=0;
  for (int n=0; n<n1-n2; n++) {
    if (b2[j] < b3[i]) {
      b1[n] = b2[j++];
    } else {
      b1[n] = b3[i++];
    }
  }
}

void psort(int n, data_t* data) {
  if (n < 32) {
    std::sort(data, data + n);
    return;
  }

  print_actual_data(data, n);

  int p = omp_get_max_threads();
  int block_size = n / p;
  int last_block_size = block_size + (n - block_size * p);
  data_t* blocks[p];

  for (int i=0; i<p; i++) {
    blocks[i] = &data[i*block_size];
  }

  printf("INITIAL BLOCK:\n\n");
  print_blocks(blocks, block_size, p, n);

#pragma omp parallel num_threads(p)
{
  int tn = omp_get_thread_num();

  int bs = (tn < p-1) ? block_size : last_block_size;
  //  printf("On thread: %d, bs=%d\n", tn, bs);
  std::sort(blocks[tn], blocks[tn] + bs);
  #pragma omp barrier
}

#pragma omp parallel for
 for (int i=0; i<p; i+=2) {
   int block_size2 = (i+1 == p-1) ? last_block_size : block_size;

   merge_blocks(blocks[i], blocks[i+1], block_size, block_size2);
 }


  // TODO: handle the left overs if max threads % 2 != 0

  printf("SORTED BLOCKS:\n\n");
  print_blocks(blocks, block_size, p, n);

  // FIXME: make this more efficient with parallelism
  std::sort(data, data + n);
}
