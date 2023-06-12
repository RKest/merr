#pragma once

#include "common.hpp"

#include "units/isq/si/mass.h"
#include "units/isq/si/thermodynamic_temperature.h"
#include "units/isq/si/energy.h"

using namespace units::isq::si::literals;

namespace r = units::isq::si::references;
namespace d = boost::math::differentiation;

// Jednostka ciepła właściwego
constexpr auto cw_ref = (r::J / (r::kg * r::K));

// Zmierzone wartości
constexpr auto g_mk = 0.180280_q_kg;
constexpr auto g_mkw = 0.372754_q_kg;
constexpr auto g_mkwl = 0.401800_q_kg;
constexpr auto g_T1 = 296.05_q_K;
constexpr auto g_T2 = 285.65_q_K;
constexpr auto T0 = 273.16_q_K;
constexpr auto Ck = 896 * cw_ref;
constexpr auto Cw = 4186 * cw_ref;

// Funkcja licząca ciepło topnienia lodu
template<typename MK, typename MKW, typename MKWL, typename T1, typename T2>
d::promote<MK, MKW, MKWL, T1, T2> qt(MK mk, MKW mkw, MKWL mkwl, T1 t1, T2 t2) {
    return ((Cw.number() * (mkw - mk) + Ck.number() * mk) * (t1 - t2) - Cw.number() * (mkwl - mkw) * (t2 - T0.number()))
           / (mkwl - mkw);
}

void lab3() {
    auto const vars = d::make_ftuple<double, 1, 1, 1, 1, 1>(g_mk.number(), g_mkw.number(),
                                                            g_mkwl.number(), g_T1.number(), g_T2.number());
    auto [mk, mkw, mkwl, t1, t2] = vars;
    auto res = qt(mk, mkw, mkwl, t1, t2);
    // Liczenie wartości ciepła topnienia lodu
    auto result = res.derivative(0, 0, 0, 0, 0) * r::J;
    // Liczehie niepewności całkowitej
    auto precision = c(
            res,
            b(0.000002),
            b(0.000002),
            b(0.000002),
            b(0.1),
            b(0.1)) * r::J;
    std::cout << fmt::format("Result:\t = {:%.2Q %q}\n", result);
    std::cout << fmt::format("Uc(qt)\t = {:%.2Q %q}\n", precision);
    std::cout << fmt::format("U\t = {:%.2Q %q}\n", precision * 2);
}