#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <random>
#include "error.h"
#include "timer.h"
#include "cp.h"

static void benchmark(std::ostream& f, int ny, int nx) {
    std::mt19937 rng;
    std::uniform_real_distribution<float> u(0.0f, 1.0f);
    std::vector<float> data(ny * nx);
    std::vector<float> result(ny * ny);
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            float v = u(rng);
            data[x + nx * y] = v;
        }
    }
    f << "cp\t" << ny << "\t" << nx << "\t" << std::flush;
    { Timer t(false, f); correlate(ny, nx, data.data(), result.data()); }
    f << std::endl;
}

static void benchmark(std::ostream& f, int ny, int nx, int iter) {
    for (int i = 0; i < iter; ++i) {
        benchmark(f, ny, nx);
    }
}

int main(int argc, const char** argv) {
    if (argc < 3 || argc > 4) {
        error("usage: cp-benchmark Y X [ITERATIONS]");
    }
    int ny = std::stoi(argv[1]);
    int nx = std::stoi(argv[2]);
    int iter = argc == 4 ? std::stoi(argv[3]) : 1;

    const char* outfile = std::getenv("PPC_OUTPUT");
    if (outfile) {
        std::ofstream f(outfile);
        if (f.fail()) {
            error(outfile, "cannot open for writing");
        }
        benchmark(f, ny, nx, iter);
        f.close();
        if (f.fail()) {
            error(outfile, "write error");
        }
    } else {
        benchmark(std::cout, ny, nx, iter);
    }
}
