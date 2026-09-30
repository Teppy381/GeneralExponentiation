// ℤ and ℝ: traits, overflow-checked integer powers, real exponents.
#pragma once

#include "power.hpp"

#include <cmath>
#include <cstdint>
#include <expected>
#include <functional>
#include <limits>
#include <system_error>
#include <type_traits>

namespace gp {

// -- Traits --------------------------------

template <Arithmetic T> struct identity_element<T, std::multiplies<>> {
    static constexpr T of(const T &) noexcept { return T(1); }
};

// Only reals have multiplicative inverses: integers stay a monoid.
template <Real T> struct inverse_operation<T, std::multiplies<>> {
    static constexpr T of(const T &x) noexcept { return T(1) / x; }
};

// Under addition every arithmetic type is a group, and power(x, n, plus) = n·x.
template <Arithmetic T> struct identity_element<T, std::plus<>> {
    static constexpr T of(const T &) noexcept { return T(0); }
};

template <Arithmetic T> struct inverse_operation<T, std::plus<>> {
    static constexpr T of(const T &x) noexcept { return static_cast<T>(-x); }
};

// -- checked_pow: the same power(), a different monoid -----------

// Multiplication lifted into std::expected: an overflow absorbs everything after it.
template <Integer T> struct checked_multiplies {
    using value_type = std::expected<T, std::errc>;

    static constexpr value_type operator()(const value_type &a,
                                           const value_type &b) noexcept {
        if (!a)
            return a;
        if (!b)
            return b;
        T r{};
        if (__builtin_mul_overflow(*a, *b, &r))
            return std::unexpected(std::errc::result_out_of_range);
        return r;
    }
};

template <Integer T>
struct identity_element<std::expected<T, std::errc>, checked_multiplies<T>> {
    static constexpr std::expected<T, std::errc>
    of(const std::expected<T, std::errc> &) noexcept {
        return T(1);
    }
};

// xⁿ with overflow reported as a value instead of UB. The check is exact:
// power() never computes a value larger than the result.
template <Integer T, Integer N>
constexpr std::expected<T, std::errc> checked_pow(T x, N n) {
    return gp::power(std::expected<T, std::errc>(x), n, checked_multiplies<T>{});
}

// -- Real exponents, built on the integer path ---------------

namespace detail {

// Integers are promoted to double, as std::pow does.
template <Arithmetic T> using promote_t = std::conditional_t<Integer<T>, double, T>;

// xʸ = xⁿ · e^{f·ln x}, where n = trunc(y) and f = y − n.
// n and f have the sign of y, so xⁿ and x^f pull the same way: xⁿ can only
// overflow or underflow when the result does. Special values follow std::pow.
template <Real R> R real_pow(R x, R y) {
    if (std::isnan(y))
        return x == R(1) ? R(1) : y;   // pow(1, NaN) == 1
    if (!(std::fabs(y) < R(0x1p63))) { // ±inf, or an integer out of int64 range
        const R ax = std::fabs(x);
        const R r = ax == R(1) ? R(1) : std::exp(y * std::log(ax));
        // Possible for long double only: a double this large is an even integer.
        const bool odd = std::isfinite(y) && std::fmod(y, R(2)) != R(0);
        return x < R(0) && odd ? -r : r;
    }
    const auto n = static_cast<std::int64_t>(y);
    const R f = y - static_cast<R>(n);
    if (f == R(0))
        return gp::pow(x, n); // integral y: the exact path, negative x included
    if (x < R(0) && std::isfinite(x))
        return std::numeric_limits<R>::quiet_NaN(); // (−8)^⅓ is not real
    const R ax = std::fabs(x); // −0 and −inf behave as their magnitudes here
    return gp::pow(ax, n) * std::exp(f * std::log(ax));
}

} // namespace detail

template <Arithmetic T, Real U> auto pow(T x, U y) {
    using R = std::common_type_t<detail::promote_t<T>, U>;
    return detail::real_pow(static_cast<R>(x), static_cast<R>(y));
}

} // namespace gp
