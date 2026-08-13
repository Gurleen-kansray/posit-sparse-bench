#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// Minimal MTX loader -> dense Eigen matrix (coordinate format, symmetric or general)
Eigen::MatrixXd load_mtx_dense(const char* f) {
    std::ifstream in(f);
    std::string line;
    std::getline(in, line); // banner
    bool symmetric = line.find("symmetric") != std::string::npos;
    while (std::getline(in, line)) if (line[0] != '%') break;
    std::istringstream dims(line);
    int rows, cols, nnz;
    dims >> rows >> cols >> nnz;
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(rows, cols);
    int r, c; double v;
    for (int k = 0; k < nnz; k++) {
        in >> r >> c >> v;
        A(r-1, c-1) = v;
        if (symmetric && r != c) A(c-1, r-1) = v;
    }
    return A;
}

int main(int argc, char** argv) {
    if (argc < 3) { std::cerr << "usage: condest_dense <matrix.mtx> <matrix_name>\n"; return 1; }
    std::string mtx_path = argv[1];
    std::string matrix_name = argv[2];

    Eigen::MatrixXd A = load_mtx_dense(mtx_path.c_str());
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A);
    double smax = svd.singularValues()(0);
    double smin = svd.singularValues()(svd.singularValues().size()-1);
    double cond = smax / smin;

    std::cout << matrix_name << ",cond_A_f64," << cond << std::endl;
    return 0;
}
