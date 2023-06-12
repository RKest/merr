#pragma once

#include "common.hpp"

#include "units/isq/si/resistance.h"

namespace r = units::isq::si::references;
namespace d = boost::math::differentiation;

// Wiadome wartości rezystorów
auto g_r1 = 39 * r::R;
auto g_r2 = 157 * r::R;
auto g_r3 = 266 * r::R;
auto g_r4 = 624 * r::R;

// Doladności przyrząów
auto res_precision = 0.1 * r::R;
auto len_precision = 0.001 * r::m;

// Funkcja liczącza rezystancję metodą Wheatstone’a
auto Xi(auto l1, auto l2, auto R) -> d::promote<decltype(l1), decltype(l2), decltype(R)> {
    return R * l1 / l2;
}

// Funkcja licząca rezystancję przy opornikach ustawionych szeregowo
auto Rsz(auto x1, auto x2, auto x3, auto x4) -> d::promote<decltype(x1), decltype(x2), decltype(x3), decltype(x4)> {
    return x1 + x2 + x3 + x4;
}

// Funkcja licząca rezystancję przy opornikach ustawionych równolegle
auto Rr(auto x1, auto x2, auto x3, auto x4) -> d::promote<decltype(x1), decltype(x2), decltype(x3), decltype(x4)> {
    return 1 / (1 / x1 + 1 / x2 + 1 / x3 + 1 / x4);
}

auto single_resistor(double prev_r, double l1, double l2, std::string_view name) {
    auto [_l1, _l2, r] = d::make_ftuple<double, 1, 1, 1>(l1, l2, prev_r);
    auto res = Xi(_l1, _l2, r);
    // Liczenie wyniku finalnie obliczonej rezystancji danego opornika
    auto result = res.derivative(0, 0, 0) * r::R;
    std::cout << fmt::format("{}: Xi\t\t = {:%.2Q %q}\n", name, result);
    // Liczenie błędu całkowitego
    auto uc = c(
            res,
            b(len_precision.number()),
            b(len_precision.number()),
            b(res_precision.number())) * r::R;
    std::cout << fmt::format("{}: Uc(Xi)\t = {:%.2Q %q}\n", name, uc);
    std::cout << fmt::format("{}: U\t\t = {:%.2Q %q}\n", name, uc * 2);
    return std::tuple{result, uc};
}

void lab1() {
    auto [calc_r1, uc1] = single_resistor(10.64, 0.595, 0.405, "R1");
    auto [calc_r2, uc2] = single_resistor(163, 0.497, 0.503, "R2");
    auto [calc_r3, uc3] = single_resistor(250.87, 0.51, 0.49, "R3");
    auto [calc_r4, uc4] = single_resistor(662.79, 0.495, 0.505, "R4");

    // Liczenie wartości rezystancji szeregowej oraz równoległej korzystając ze znanych wartości rezystorów, oraz
    // przyjmując niedokładność pomiaru tych rezystorów na taką, jaka jest dokładność opornika dekadowego
    {
        auto [r1, r2, r3, r4] = d::make_ftuple<double, 1, 1, 1, 1>(g_r1.number(), g_r2.number(), g_r3.number(),
                                                                   g_r4.number());
        auto res_sz = Rsz(r1, r2, r3, r4);
        auto result_sz = res_sz.derivative(0, 0, 0, 0) * r::R;
        auto res_r = Rr(r1, r2, r3, r4);
        auto result_r = res_r.derivative(0, 0, 0, 0) * r::R;
        auto uc_sz = c(res_sz,
                       b(res_precision).number(),
                       b(res_precision).number(),
                       b(res_precision).number(),
                       b(res_precision).number()) * r::R;
        auto uc_r = c(res_r,
                      b(res_precision).number(),
                      b(res_precision).number(),
                      b(res_precision).number(),
                      b(res_precision).number()) * r::R;
        std::cout << fmt::format("From precise values: Rsz\t = {:%.2Q %q}\n", result_sz);
        std::cout << fmt::format("From precise values: Rr\t\t = {:%.2Q %q}\n", result_r);
        std::cout << fmt::format("From precise values: Uc(Rsz)\t = {:%.2Q %q}\n", uc_sz);
        std::cout << fmt::format("From precise values: U\t\t = {:%.2Q %q}\n", uc_sz * 2);
        std::cout << fmt::format("From precise values: Uc(Rr)\t = {:%.2Q %q}\n", uc_r);
        std::cout << fmt::format("From precise values: U\t\t = {:%.2Q %q}\n", uc_r * 2);
    }

    // Liczenie tego co wyżej, jedynie korzystając z wartości rezystorów obliczonych za pomocą metody Whitestone'a
    {
        auto [r1, r2, r3, r4] = d::make_ftuple<double, 1, 1, 1, 1>(calc_r1.number(), calc_r2.number(), calc_r3.number(),
                                                                   calc_r4.number());
        auto res_sz = Rsz(r1, r2, r3, r4);
        auto result_sz = res_sz.derivative(0, 0, 0, 0) * r::R;
        auto res_r = Rr(r1, r2, r3, r4);
        auto result_r = res_r.derivative(0, 0, 0, 0) * r::R;
        auto uc_sz = c(res_sz,
                       b(uc1).number(),
                       b(uc2).number(),
                       b(uc3).number(),
                       b(uc4).number()) * r::R;
        auto uc_r = c(res_r,
                      b(uc1).number(),
                      b(uc2).number(),
                      b(uc3).number(),
                      b(uc4).number()) * r::R;
        std::cout << fmt::format("From calculated values: Rsz\t = {:%.2Q %q}\n", result_sz);
        std::cout << fmt::format("From calculated values: Rr\t = {:%.2Q %q}\n", result_r);
        std::cout << fmt::format("From calculated values: Uc(Rsz)\t = {:%.2Q %q}\n", uc_sz);
        std::cout << fmt::format("From calculated values: U\t = {:%.2Q %q}\n", uc_sz * 2);
        std::cout << fmt::format("From calculated values: Uc(Rr)\t = {:%.2Q %q}\n", uc_r);
        std::cout << fmt::format("From calculated values: U\t = {:%.2Q %q}\n", uc_r * 2);
    }
}