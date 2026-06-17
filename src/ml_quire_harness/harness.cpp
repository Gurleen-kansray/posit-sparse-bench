// harness.cpp: four-way dot product comparison on REAL trained MLP weights/activations
// methods: float32 (rounded accumulation), double64 (rounded accumulation),
//          posit32 (rounded accumulation, no quire), posit32+quire (fused dot product)
#include <universal/number/posit/posit.hpp>
#include <blas/blas.hpp>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>

using namespace sw::universal;
using namespace sw::numeric::containers;
using namespace sw::blas;

int main() {
    std::ifstream fin("dotpairs.txt");
    if (!fin) { std::cerr << "cannot open dotpairs.txt\n"; return 1; }

    size_t total;
    fin >> total;

    std::ofstream fout("harness_results.csv");
    fout << std::setprecision(17);
    fout << "layer,float32,double64,posit32,posit32_quire\n";

    for (size_t p = 0; p < total; ++p) {
        std::string tag;
        size_t n;
        fin >> tag >> n;

        std::vector<double> xin(n), win(n);
        for (size_t i = 0; i < n; ++i) fin >> xin[i];
        for (size_t i = 0; i < n; ++i) fin >> win[i];

        // float32 path
        {
            vector<float> x(n), w(n);
            for (size_t i = 0; i < n; ++i) { x[i] = static_cast<float>(xin[i]); w[i] = static_cast<float>(win[i]); }
            float r = dot(x, w);
            fout << tag << "," << r << ",";
        }
        // double64 path
        {
            vector<double> x(n), w(n);
            for (size_t i = 0; i < n; ++i) { x[i] = xin[i]; w[i] = win[i]; }
            double r = dot(x, w);
            fout << r << ",";
        }
        // posit32 (no quire) path
        {
            using P32 = posit<32,2>;
            vector<P32> x(n), w(n);
            for (size_t i = 0; i < n; ++i) { x[i] = xin[i]; w[i] = win[i]; }
            P32 r = dot(x, w);
            fout << r << ",";
        }
        // posit32 + quire (fused dot product) path
        {
            using P32 = posit<32,2>;
            vector<P32> x(n), w(n);
            for (size_t i = 0; i < n; ++i) { x[i] = xin[i]; w[i] = win[i]; }
            P32 r = fdp(x, w);
            fout << r << "\n";
        }
    }
    std::cerr << "wrote harness_results.csv (" << total << " rows)\n";
    return 0;
}
