#include "cp.h"
#include <math.h>
#include <stdio.h>

double get_mean(const float* row, int nx) {
    double sum = 0.0, mean;
    for (int x=0; x<nx; x++) {
	sum += (double) row[x];
    }
    mean = sum / (double) nx;
    return mean;
}

double get_root_square_sum(const float* row, int nx) {
    double square_sum = 0, root_square_sum;
    for (int x=0; x<nx; x++) {
	square_sum += pow((double) row[x], 2.0);
    }
    root_square_sum = sqrt(square_sum);
    return root_square_sum;
}

void correlate(int ny, int nx, const float* data, float* result) {
    int dptr;
    double row_mean, row_rss;
    double X[nx*ny];

    printf("\n");
    for (int y=0; y<ny; y++) {
	row_mean = get_mean(&data[nx*y], nx);
	row_rss = get_root_square_sum(&data[nx*y], nx);
	printf("row: ");
	for (int x=0; x<nx; x++) {
	    dptr = nx*y + x;
	    X[dptr] = (data[dptr] - row_mean) / row_rss;
	    result[dptr] = X[dptr];
	    printf("%f ", X[dptr]);
	}
	printf("row mean: %f, row rss: %f", row_mean, row_rss);
	printf("\n");
    }
}
