#pragma once

#include <numeric>
#include <cmath>
#include <array>
#include <iostream>

#include <units/isq/si/length.h>
#include <boost/math/differentiation/autodiff.hpp>
#include <units/format.h>

auto a(auto iterable) -> decltype(iterable)::value_type {
    using value_t = typename decltype(iterable)::value_type;
    auto const n = static_cast<double>(iterable.size());
    auto const avg = std::accumulate(begin(iterable), end(iterable), value_t{}) / iterable.size();
    auto const acc = std::accumulate(begin(iterable), end(iterable), value_t{},
                                     [&avg](auto a, auto b) {
                                         return a + value_t{(b.number() - avg.number()) * (b.number() - avg.number())};
                                     });
    return value_t{std::sqrt((1 / n * (n - 1)) * acc.number())};
}

auto b(auto single) -> decltype(single)
{
    static auto root_3 = sqrt(3);
    return single / root_3;
}
