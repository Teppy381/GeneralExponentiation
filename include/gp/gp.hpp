// The umbrella header. Include this one only: a trait specialization has to be
// visible everywhere a concept is checked, or programs become ill-formed (NDR).
#pragma once

#include "concepts.hpp" // layer 0: concepts and customization points
#include "power.hpp"    // layer 1: power(x, n, op) and pow(x, n)
#include "scalar.hpp"   // ℤ, ℝ: traits, checked_pow, real exponents
#include "complex.hpp"  // ℂ: complex exponents
#include "matrix.hpp"   // Eigen: Aⁿ and Aᵗ
