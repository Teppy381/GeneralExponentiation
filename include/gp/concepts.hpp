// Layer 0: what a type has to provide to be raised to a power.
#pragma once

#include <concepts>
#include <functional>
#include <string_view>
#include <type_traits>

namespace gp {

// -- Domain concepts -----------------------------

template <class N>
concept Integer = std::integral<N> && !std::same_as<std::remove_cv_t<N>, bool>;

template <class T>
concept Real = std::floating_point<T>;

template <class T>
concept Arithmetic = Integer<T> || Real<T>;

// Structural: std::complex and anything shaped like it.
template <class C>
concept Complex = Real<typename C::value_type> && requires(const C &z) {
    { z.real() } -> std::same_as<typename C::value_type>;
    { z.imag() } -> std::same_as<typename C::value_type>;
};

// Every non-zero element has a multiplicative inverse.
template <class S>
concept Field = Real<S> || Complex<S>;

// -- Customization points --------------------------
// Specialize them for (T, Op) to adapt a type without touching it.

// Provides `static T of(const T& sample)`, the identity element of Op.
// It takes a sample because the identity may depend on a run-time shape (MatrixXd).
template <class T, class Op> struct identity_element {};

// Provides `static T of(const T& x)`, the inverse of x with respect to Op.
template <class T, class Op> struct inverse_operation {};

// For types that are closed under Op syntactically but not mathematically
// (non-square matrices). The same idiom as std::ranges::disable_sized_range.
template <class T, class Op> inline constexpr bool disable_semigroup = false;

// -- Algebraic concepts --------------------------
// Each concept refines the previous one, so overloads constrained by them are
// ordered by subsumption. The axioms are semantic: the compiler cannot check them.

// Axiom: op is associative.
template <class T, class Op>
concept Semigroup =
    !disable_semigroup<T, Op> && std::copyable<T> &&
    std::regular_invocable<Op &, const T &, const T &> &&
    // convertible_to, not same_as: Eigen returns expression templates.
    std::convertible_to<std::invoke_result_t<Op &, const T &, const T &>, T>;

// Axiom: op(e, x) == op(x, e) == x.
template <class T, class Op>
concept Monoid = Semigroup<T, Op> && requires(const T &x) {
    { identity_element<T, Op>::of(x) } -> std::convertible_to<T>;
};

// Axiom: op(x, x⁻¹) == op(x⁻¹, x) == e.
template <class T, class Op>
concept Group = Monoid<T, Op> && requires(const T &x) {
    { inverse_operation<T, Op>::of(x) } -> std::convertible_to<T>;
};

template <class T, class Op>
    requires Monoid<T, Op>
constexpr T identity(const T &sample, const Op &) {
    return identity_element<T, Op>::of(sample);
}

template <class T, class Op>
    requires Group<T, Op>
constexpr T inverse(const T &x, const Op &) {
    return inverse_operation<T, Op>::of(x);
}

// The strongest structure (T, Op) models.
template <class T, class Op = std::multiplies<>>
consteval std::string_view structure_of() {
    if constexpr (Group<T, Op>)
        return "Group";
    else if constexpr (Monoid<T, Op>)
        return "Monoid";
    else if constexpr (Semigroup<T, Op>)
        return "Semigroup";
    else
        return "none";
}

} // namespace gp
