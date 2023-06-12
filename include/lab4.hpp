#pragma once

#include "common.hpp"

#include <algorithm>

#include "units/isq/si/length.h"
#include "units/isq/si/acceleration.h"
#include "units/isq/si/mass.h"
#include "fmt/format.h"

using namespace units::isq::si::literals;

// Stałe oraz wartości zmierzone
constexpr auto g_const = 9.81_q_m_per_s2;
constexpr auto x0_1 = 0.16_q_m;
constexpr auto x0_2 = 0.243_q_m;
constexpr std::array weights{50.552_q_g, 50.396_q_g, 50.682_q_g, 50.278_q_g, 50.402_q_g, 50.330_q_g, 50.554_q_g};
constexpr std::array xs_1{0.195_q_m, 0.228_q_m, 0.261_q_m, 0.295_q_m, 0.328_q_m, 0.36_q_m, 0.394_q_m};
constexpr std::array xs_2{0.264_q_m, 0.285_q_m, 0.307_q_m, 0.329_q_m, 0.351_q_m, 0.372_q_m, 0.393_q_m};

// Funkcja licząca stałą sprężyny za pomocą prawa Hooke'a
template<typename M, typename X>
auto k(const M &m, const X &x) {
    return m * g_const / x;
}

// Sumowanie pierwszych 'n' ciężarków
auto sum_of_weights_until(int up_unit) {
    return std::accumulate(begin(weights), begin(weights) + up_unit + 1, 0.0_q_g);
}

// Liczenie dla wszystkich pomiarów stałej sprężyny
auto calc(decltype(xs_1) xs, decltype(x0_1) x0) {
    std::vector<decltype(k(1_q_g, 1_q_m))> ress;
    for (int i = 0; i < weights.size(); ++i) {
        auto const m = sum_of_weights_until(i);
        auto const x = xs.at(i) - x0;
        auto const res = k(m, x);
        ress.push_back(res);
    }
    return ress;
}

void lab4() {
    auto res1 = calc(xs_1, x0_1);
    auto res2 = calc(xs_2, x0_2);
    // Liczenie niepewności całkowitej jako niepewność typu A
    auto err1 = a(res1);
    auto err2 = a(res2);
    std::cout << fmt::format("Uc(k) for spring 1\t = {:%.2Q %q}\n", err1);
    std::cout << fmt::format("Uc(k) for spring 2\t = {:%.2Q %q}\n", err2);
    std::cout << fmt::format("U for spring 1\t\t = {:%.2Q %q}\n", err1 * 2);
    std::cout << fmt::format("U for spring 2\t\t = {:%.2Q %q}\n", err2 * 2);
}
