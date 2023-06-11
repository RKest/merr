#pragma once

#include "common.hpp"

#include "units/isq/si/mass.h"
#include "units/isq/si/thermodynamic_temperature.h"
#include "units/isq/si/energy.h"

using namespace units::isq::si::literals;

namespace r = units::isq::si::references;
namespace d = boost::math::differentiation;

constexpr auto cw_ref = (r::J / (r::kg * r::K));

constexpr auto g_mk = 0.180280_q_kg;
constexpr auto g_mkw = 0.372754_q_kg;
constexpr auto g_mkwl = 0.401800_q_kg;
constexpr auto g_T1 = 296.05_q_K;
constexpr auto g_T2 = 285.65_q_K;
constexpr auto T0 = 273.16_q_K;
constexpr auto Ck = 896 * cw_ref;
constexpr auto Cw = 4186 * cw_ref;

template<typename MK, typename MKW, typename MKWL, typename T1, typename T2>
d::promote<MK, MKW, MKWL, T1, T2> qt(MK mk, MKW mkw, MKWL mkwl, T1 t1, T2 t2) {
    return ((Cw.number() * (mkw - mk) + Ck.number() * mk) * (t1 - t2) - Cw.number() * (mkwl - mkw) * (t2 - T0.number()))
           / (mkwl - mkw);
}

void lab3() {
    constexpr unsigned max_order = 1;
    auto const vars = d::make_ftuple<double, max_order, max_order, max_order, max_order, max_order>(g_mk.number(),
                                                                                                    g_mkw.number(),
                                                                                                    g_mkwl.number(),
                                                                                                    g_T1.number(),
                                                                                                    g_T2.number());
    auto [_mk, _mkw, _mkwl, _t1, _t2] = vars;
    auto res = qt(_mk, _mkw, _mkwl, _t1, _t2);
    std::cout << std::setprecision(20) << "Res: " << res.derivative(0, 0, 0, 0, 0) << '\n';
}