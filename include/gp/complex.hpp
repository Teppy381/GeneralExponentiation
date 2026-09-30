// ℂ: traits, real and complex exponents, built on the real layer.
#pragma once

#include "scalar.hpp"

#include <cmath>
#include <complex>
#include <cstdint>
#include <functional>
#include <limits>

namespace gp {

// -- Traits: ℂ under multiplication is a group (0 aside) -----------

template <Complex C> struct identity_element<C, std::multiplies<>> {
    static constexpr C of(const C &) {
        return C(typename C::value_type(1), typename C::value_type(0));
    }
};

template <Complex C> struct inverse_operation<C, std::multiplies<>> {
    static constexpr C of(const C &z) {
        return identity_element<C, std::multiplies<>>::of(z) / z;
    }
};

// -- Real and complex exponents ----------------------

namespace detail {

template <Complex C> using value_t = typename C::value_type;

// zʷ = exp(w · Log z) for w = a + ib on the principal branch:
// |z|ᵃ · e^{−bθ} · cis(b·ln|z| + aθ), θ = arg z. |z|ᵃ comes from the real layer.
template <Complex C> C complex_pow(const C &z, value_t<C> a, value_t<C> b) {
    using R = value_t<C>;
    const R r = std::hypot(z.real(), z.imag());
    if (r == R(0)) { // 0ʷ is 0 for Re w > 0 and undefined otherwise (w ≠ 0 here)
        const R nan = std::numeric_limits<R>::quiet_NaN();
        return a > R(0) ? C(R(0), R(0)) : C(nan, nan);
    }
    const R theta = std::atan2(z.imag(), z.real());
    const R rho = gp::pow(r, a) * std::exp(-b * theta);
    const R phi = b * std::log(r) + a * theta;
    return C(rho * std::cos(phi), rho * std::sin(phi));
}

} // namespace detail

// An integral y takes the exact route through power(): (1+i)⁸ == 16.
template <Complex C, Real U> C pow(const C &z, U y) {
    using R = detail::value_t<C>;
    const auto e = static_cast<R>(y);
    if (std::trunc(e) == e && std::fabs(e) < R(0x1p63))
        return gp::pow(z, static_cast<std::int64_t>(e));
    return detail::complex_pow(z, e, R(0));
}

template <Complex C, Complex W> C pow(const C &z, const W &w) {
    using R = detail::value_t<C>;
    if (w.imag() == 0)
        return gp::pow(z, static_cast<R>(w.real()));
    return detail::complex_pow(z, static_cast<R>(w.real()), static_cast<R>(w.imag()));
}

template <Arithmetic T, Complex W> W pow(T x, const W &w) {
    using R = detail::value_t<W>;
    return gp::pow(W(static_cast<R>(x), R(0)), w);
}

} // namespace gp
