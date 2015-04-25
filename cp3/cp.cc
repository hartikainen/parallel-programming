#include "cp.h"
#include <cmath>
#include <stdio.h>
#include <iostream>
#include "../common/vector.h"

double get_mean(const float* row, int nx) {
    double sum = 0.0, mean;
    for (int x=0; x<nx; x++) {
	sum += (double) row[x];
    }
    mean = sum / (double) nx;
    return mean;
}

double get_root_square_sum(double* row, int nx) {
    double square_sum = 0, root_square_sum;
    for (int x=0; x<nx; x++) {
	square_sum += pow((double) row[x], 2.0);
    }
    root_square_sum = std::sqrt(square_sum);
    return root_square_sum;
}

double get_row_sum(double4_t* row, int n) {
    double4_t sum4_t;
    double sum = 0.0;
    for (int i=0; i<n; i++) {
	sum4_t += 0.0;
    }

    return sum;
}

void correlate(int ny, int nx, const float* data, float* result) {
    int nnx = std::ceil(nx/4.0);
    int xpad = nnx * 4.0 - nx;
    double4_t* X = double4_alloc(nnx*ny);

    // Copy the data into double4_t array X
    for (int y=0; y<ny; y++) {
	for (int x1=0; x1<nnx-1; x1++) {
	    for (int x2=0; x2<4; x2++) {
		X[y*nnx + x1][x2] = data[y*nx + x1*4 + x2];
		printf("times here: %d", x2);
	    }
	}

	for (int x=0; x<xpad; x++) {
	    X[y*nnx + nnx-1][x] = data[y*nx + (nnx - 1) * 4 + x];
	}
	// Initialize the remaining 'right' padding to 0.0
	for (int x=4-xpad; x<4; x++) {
	    X[y*nnx + nnx-1][x] = 0.0;
	}
    }

    std::cout << "\ndata: \n";
    for (int y=0; y<ny; y++) {
	for (int x=0; x<nx; x++) {
	    printf("%f ", data[y*nx + x]);
	}
    	printf("\n");
    }

    std::cout << "\nX: \n";
    for (int y=0; y<ny; y++) {
	for (int x1=0; x1<nnx; x1++) {
	    for (int x2=0; x2<4; x2++) {
		printf("%f ", X[y*nnx + x1][x2]);
	    }
	}
	printf("\n");
    }

// #pragma omp parallel for
//     for (int y=0; y<ny; y++) {
// 	double row_mean = get_mean(&data[y*nx], nx);

// 	for (int x=0; x<nx; x++) {
// 	    X[y*nx + x] = ((double) data[y*nx + x]) - row_mean;

// 	}

// 	double row_rss = get_root_square_sum(&X[y*nx], nx);

// 	for (int x=0; x<nx; x++) {
// 	    X[y*nx + x] = X[y*nx + x] / row_rss;
// 	}
//     }

// #pragma omp parallel for schedule(static, 1)
//     for (int y=0; y<ny; y++) {
// 	for (int x=y; x<ny; x++) {
// 	    double r = 0.0;
// 	    for (int i=0; i<nx; i++) {
// 		r += X[x*nx + i] * X[y*nx + i];
// 	    }
// 	    result[y*ny + x] = r;
// 	}
//     }
}
