#include "cp.h"
#include <cuda_runtime.h>

void correlate(int ny, int nx, const float* data, float* result) {
    // FIXME
    for (int i = 0; i < ny * ny; ++i) {
        result[i] = 0.0f;
    }
}
