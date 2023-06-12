#pragma once

#include "common.hpp"

using namespace units::isq::si::literals;

// Wartości długości ogniskowych
std::array fs{
        0.082648115_q_m,
        0.083753894_q_m,
        0.084975719_q_m,
        0.082099663_q_m,
        0.082368789_q_m,
        0.083969321_q_m,
        0.082580077_q_m,
        0.083305845_q_m,
        0.085184719_q_m,
        0.084529543_q_m,
};

void lab2() {
    // Liczenie niepewności całkowitej jako niepewność typu A
    const auto res = a(fs);
    std::cout << fmt::format("Uc(f)\t = {:%.2GQ %q}\n", res);
    std::cout << fmt::format("U\t = {:%.2GQ %q}\n", res * 2);
}