#include "mf.h"
#include <cuda_runtime.h>

void mf(int ny, int nx, int hy, int hx, const float* in, float* out) {
    // FIXME
    for (int i = 0; i < ny * nx; ++i) {
        out[i] = in[i];
    }
}
