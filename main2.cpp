#include <iostream>
#include <algorithm>
#include <numeric>

#include "boost/math/differentiation/autodiff.hpp"
#include "boost/multiprecision/cpp_bin_float.hpp"

namespace d = boost::math::differentiation;

using float50 = boost::multiprecision::cpp_bin_float_50;
constexpr unsigned max_order = 1;
constexpr double g = 9.81;
constexpr double x0_m = 0.16;
constexpr double x02_m = 0.243;
constexpr std::array weights{50.552, 50.396, 50.682, 50.278, 50.402, 50.330, 50.554};
constexpr std::array xs_m{0.195, 0.228, 0.261, 0.295, 0.328, 0.36, 0.394};
constexpr std::array xs2_m{0.264, 0.285, 0.307, 0.329, 0.351, 0.372, 0.393};

template<typename M, typename X>
d::promote<M, X> k(const M &m, const X &x) {
    return m * g / x;
}

auto sum_of_weights_until(int up_unit) {
    return std::accumulate(begin(weights), begin(weights) + up_unit + 1, 0.0);
}

void calc(decltype(xs_m) local_xs_m, double x0) {
    for (int i = 0; i < weights.size(); ++i) {
        auto const summed_weights = sum_of_weights_until(i);
        auto const vars = d::make_ftuple<float50, max_order, max_order>(summed_weights / 1000.0, local_xs_m.at(i) - x0);
        auto const m = std::get<0>(vars);
        auto const x = std::get<1>(vars);
        auto const res = k(m, x);
        std::cout << "Res = " << res.derivative(0, 0) << '\t';
        std::cout << "Derivative over m = " << res.derivative(1, 0) << '\t';
        std::cout << "Derivative over x = " << res.derivative(0, 1) << '\n';
    }
}

auto average(const auto &iterable) -> decltype(iterable)::element_type {
    return std::accumulate(begin(iterable), end(iterable)) / iterable.size();
}

auto a(const auto &iterable) -> decltype(iterable)::element_type {
    auto const n = iterable.size();
    auto const avg = average(iterable);
    auto const acc = std::accumulate(begin(iterable), end(iterable), decltype(iterable)::element_type(),
                                     [&avg](auto a, auto b) {
                                         return a + std::pow(b - avg, 2);
                                     });
    return 1 / n * (n - 1) * acc;
}


int main() {
    calc(xs_m, x0_m);
}
