// Layer 1: Stepanov's power algorithm, one overload per algebraic structure.
#pragma once

#include "concepts.hpp"

#include <cstdio>
#include <cstdlib>
#include <functional>
#include <source_location>
#include <type_traits>
#include <utility>

namespace gp::detail {

// Aborts at run time. During constant evaluation a call to a non-constexpr
// function is ill-formed, so a violated precondition becomes a compile error.
[[noreturn]] inline void precondition_failed(const char* condition,
                                             std::source_location where = std::source_location::current()) {
    std::fprintf(stderr, "%s:%u: %s: precondition failed: %s\n", where.file_name(),
                 static_cast<unsigned>(where.line()), where.function_name(), condition);
    std::abort();
}

}  // namespace gp::detail

#define GP_EXPECTS(...) ((__VA_ARGS__) ? void(0) : ::gp::detail::precondition_failed(#__VA_ARGS__))

namespace gp {

namespace detail {

template <Integer N>
constexpr bool odd(N n) noexcept {
    return (n & N(1)) != N(0);
}

template <Integer N>
constexpr N half(N n) noexcept {  // n >= 0
    return static_cast<N>(n >> 1);
}

template <Integer N>
constexpr bool is_negative(N n) noexcept {
    if constexpr (std::is_signed_v<N>)
        return n < N(0);
    else
        return false;
}

// |n| as an unsigned number: well-defined even for the most negative n.
template <Integer N>
constexpr auto magnitude(N n) noexcept {
    using U = std::make_unsigned_t<N>;
    return is_negative(n) ? static_cast<U>(U(0) - static_cast<U>(n)) : static_cast<U>(n);
}

// Calls op exactly the way the concepts checked it and materializes the result
// as T: this evaluates Eigen expression templates and undoes integer promotion.
template <class T, class Op>
constexpr T apply(Op& op, const T& a, const T& b) {
    return static_cast<T>(op(a, b));
}

// r·xⁿ for n >= 0. Invariant: r·xⁿ is the same on every iteration.
// Returns right after the last multiplication: no squaring past the top bit,
// so no intermediate value is larger than the result.
template <class T, Integer N, class Op>
constexpr T power_accumulate(T r, T x, N n, Op& op) {
    if (n == N(0)) return r;
    while (true) {
        if (odd(n)) {
            r = apply(op, r, x);
            if (n == N(1)) return r;
        }
        n = half(n);
        x = apply(op, x, x);
    }
}

// xⁿ for n > 0: no identity element needed.
template <class T, Integer N, class Op>
constexpr T power_semigroup(T x, N n, Op& op) {
    while (!odd(n)) {
        x = apply(op, x, x);
        n = half(n);
    }
    if (n == N(1)) return x;
    return power_accumulate(x, apply(op, x, x), half(static_cast<N>(n - N(1))), op);
}

}  // namespace detail

// ---- power(x, n, op): Semigroup ⊂ Monoid ⊂ Group, chosen by subsumption ----

template <class T, Integer N, class Op = std::multiplies<>>
    requires Semigroup<T, Op>
constexpr T power(T x, N n, Op op = {}) {
    GP_EXPECTS(n > N(0));  // a semigroup has no x⁰
    return detail::power_semigroup(std::move(x), n, op);
}

template <class T, Integer N, class Op = std::multiplies<>>
    requires Monoid<T, Op>
constexpr T power(T x, N n, Op op = {}) {
    GP_EXPECTS(!detail::is_negative(n));  // a monoid has no inverses
    if (n == N(0)) return gp::identity(x, op);
    return detail::power_semigroup(std::move(x), n, op);
}

template <class T, Integer N, class Op = std::multiplies<>>
    requires Group<T, Op>
constexpr T power(T x, N n, Op op = {}) {
    if (n == N(0)) return gp::identity(x, op);
    if (detail::is_negative(n)) return detail::power_semigroup(gp::inverse(x, op), detail::magnitude(n), op);
    return detail::power_semigroup(std::move(x), n, op);
}

// ---- pow(x, n): any multiplicative semigroup, integer exponent ------------
// One overload serves integers, reals, complex numbers and matrices: the traits
// of T decide whether n may be zero or negative.
template <class T, Integer N>
    requires Semigroup<T, std::multiplies<>>
constexpr T pow(const T& x, N n) {
    return gp::power(x, n, std::multiplies<>{});
}

}  // namespace gp
