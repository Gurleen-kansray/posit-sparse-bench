#include <universal/number/posit/posit.hpp>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
using namespace sw::universal;
int main(int argc, char** argv) {
  for (int a = 1; a < argc; ++a) {
    std::ifstream in(argv[a]); std::string line;
    do { std::getline(in, line); } while (!line.empty() && line[0] == '%');
    long r, c, nz; std::istringstream(line) >> r >> c >> nz;
    long exact = 0, total = 0; double maxrel = 0;
    while (std::getline(in, line)) {
      std::istringstream s(line); long i, j; double v;
      if (!(s >> i >> j >> v)) continue;
      posit<32,2> p(v); double back = double(p);
      ++total; if (back == v) ++exact;
      if (v != 0) maxrel = std::max(maxrel, std::fabs(back - v) / std::fabs(v));
    }
    std::cout << a << " " << argv[a] << ": " << exact << "/" << total
              << " exact, max rel err " << maxrel << "\n";
  }
}
