#pragma once

#include <numeric>
#include <cmath>
#include <array>
#include <iostream>

#include <units/isq/si/length.h>
#include <boost/math/differentiation/autodiff.hpp>
#include <units/format.h>

// Funkcja licząca niepewność typu A
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

// Funkcja licząca niepewność typu B
auto b(auto single) -> decltype(single) {
    static auto root_3 = sqrt(3);
    return single / root_3;
}

// Detale implementacyjne
template<typename ...Ts>
struct type_list {
};

template<std::size_t Sz, std::size_t I>
auto make_mask() {
    return []<std::size_t ...Is>(std::index_sequence<Is...>) {
        return std::index_sequence<static_cast<std::size_t>(Is == I)...>{};
    }(std::make_index_sequence<Sz>());
}

template<std::size_t Sz>
auto make_all_masks() {
    return []<std::size_t ...Is>(std::index_sequence<Is...>) {
        return type_list<decltype(make_mask<Sz, Is>())...>();
    }(std::make_index_sequence<Sz>());
}

// Funkcja licząca jeden kawałek niepewności całkowitej
template<std::size_t ...Is>
auto c_impl(auto res, auto u, std::index_sequence<Is...>) {
    return pow(res.derivative(Is...) * u, 2);
}

// Funkcja licząca niepewność całkowitą
auto c(auto autodiff_result, auto... measurement_uncertainties) {
    return std::sqrt([&]<typename ...Seq>(type_list<Seq...>) {
        return (... + c_impl(autodiff_result, measurement_uncertainties, Seq{}));
    }(make_all_masks<sizeof...(measurement_uncertainties)>()));
}