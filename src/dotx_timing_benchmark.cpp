#include <cstdio>
#include <cstdint>
#include <vector>
#include <chrono>
#include <random>
#include <utility>
#include <universal/number/posit/posit.hpp>
#include <universal/number/posit/fdp.hpp>
using namespace sw::universal;

double dot_double(const std::vector<double>& a, const std::vector<double>& b){
    double s = 0;
    for (size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
    return s;
}

double dot_posit32_naive(const std::vector<double>& a, const std::vector<double>& b){
    posit<32,2> s = 0;
    for (size_t i = 0; i < a.size(); ++i){
        posit<32,2> pa(a[i]), pb(b[i]);
        s = s + pa * pb;
    }
    return double(s);
}

double dot_posit32_quire(const std::vector<double>& a, const std::vector<double>& b){
    quire<32,2,2> q = 0;
    for (size_t i = 0; i < a.size(); ++i){
        posit<32,2> pa(a[i]), pb(b[i]);
        q += quire_mul(pa, pb);
    }
    posit<32,2> pr; convert(q.to_value(), pr);
    return double(pr);
}

// Time-bounded: run each function for a fixed wall-clock budget and count
// how many calls fit, rather than guessing a call count up front. Keeps
// total runtime predictable regardless of how slow posit/quire turn out
// to be relative to native double arithmetic.
template<typename F>
std::pair<double,long> time_ns_per_call(F&& f, const std::vector<double>& a, const std::vector<double>& b, double budget_seconds){
    volatile double sink = 0;
    for (int r = 0; r < 5; ++r){ sink += f(a,b); asm volatile("" : : : "memory"); }

    long count = 0;
    auto t0 = std::chrono::high_resolution_clock::now();
    auto deadline = t0 + std::chrono::duration<double>(budget_seconds);
    while (std::chrono::high_resolution_clock::now() < deadline){
        sink += f(a,b);
        asm volatile("" : : : "memory");
        ++count;
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    (void)sink;
    double ns_per_call = std::chrono::duration<double,std::nano>(t1 - t0).count() / (double)count;
    return {ns_per_call, count};
}

int main(){
    struct Case { int n; const char* label; };
    std::vector<Case> cases = {
        {112,  "bcsstk03 (112)"},
        {1806, "bcsstk14 (1806)"},
        {4960, "add32 (4960)"}
    };
    const double BUDGET = 1.0; // seconds per function per case

    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    printf("%-18s %-12s %-16s %-16s %-10s %-10s %-12s %-12s\n",
           "matrix(n)","double_ns","posit32naive_ns","posit32quire_ns",
           "n_dbl","n_quire","quire/dbl","quire/naive");

    for (auto& c : cases){
        std::vector<double> a(c.n), b(c.n);
        for (int i = 0; i < c.n; ++i){ a[i] = dist(rng); b[i] = dist(rng); }

        auto [t_double, n_double] = time_ns_per_call(dot_double, a, b, BUDGET);
        auto [t_naive,  n_naive]  = time_ns_per_call(dot_posit32_naive, a, b, BUDGET);
        auto [t_quire,  n_quire]  = time_ns_per_call(dot_posit32_quire, a, b, BUDGET);

        printf("%-18s %-12.2f %-16.2f %-16.2f %-10ld %-10ld %-12.2f %-12.2f\n",
               c.label, t_double, t_naive, t_quire,
               n_double, n_quire,
               t_quire / t_double, t_quire / t_naive);
    }
    return 0;
}
