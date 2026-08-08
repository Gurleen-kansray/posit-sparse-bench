#include <universal/number/posit/posit.hpp>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <iostream>
#include <fstream>
#include <string>

using namespace sw::universal;

template<size_t nbits>
Eigen::SparseMatrix<double> quantize(const Eigen::SparseMatrix<double>& A) {
    Eigen::SparseMatrix<double> Aq = A;
    for (int k = 0; k < Aq.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(Aq, k); it; ++it) {
            posit<nbits, 2> p = it.value();
            it.valueRef() = double(p);
        }
    return Aq;
}

double condest_f64(const Eigen::SparseMatrix<double>& A) {
    // reuse your existing CMSW/condest float64 implementation here
    // return condition number estimate
}

int main(int argc, char** argv) {
    std::string mtx_path = argv[1];
    std::string matrix_name = argv[2];

    Eigen::SparseMatrix<double> A = loadMarketMatrix(mtx_path); // match your existing loader name

    for (int bits : {8, 16, 32, 64}) {
        Eigen::SparseMatrix<double> Aq;
        if (bits == 8)  Aq = quantize<8>(A);
        if (bits == 16) Aq = quantize<16>(A);
        if (bits == 32) Aq = quantize<32>(A);
        if (bits == 64) Aq = quantize<64>(A);

        double c = condest_f64(Aq);
        std::cout << matrix_name << "," << bits << "," << c << std::endl;
    }
    return 0;
}
