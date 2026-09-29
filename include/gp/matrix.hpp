// Eigen matrices: integer powers through the traits, real and complex exponents
// through the spectral decomposition, built on the complex layer.
#pragma once

#include "complex.hpp"

#include <cmath>
#include <complex>
#include <concepts>
#include <functional>
#include <limits>
#include <stdexcept>

#include <Eigen/Dense>

namespace gp {

// ---- Concepts --------------------------------------------------------------

// A plain Eigen matrix: not an expression template and not an Array.
template <class M>
concept EigenMatrix = requires {
    typename M::Scalar;
    typename M::PlainObject;
} && std::derived_from<M, Eigen::PlainObjectBase<M>> && std::derived_from<M, Eigen::MatrixBase<M>>;

// Square at compile time, or sized at run time (then checked by assertions).
template <class M>
concept SquareMatrix = EigenMatrix<M> && (M::RowsAtCompileTime == M::ColsAtCompileTime ||
                                          M::RowsAtCompileTime == Eigen::Dynamic ||
                                          M::ColsAtCompileTime == Eigen::Dynamic);

// ---- Traits ----------------------------------------------------------------

// Square matrices over a ring form a monoid; e is the identity of the sample's size.
template <SquareMatrix M>
struct identity_element<M, std::multiplies<>> {
    static M of(const M& x) {
        GP_EXPECTS(x.rows() == x.cols());
        return M::Identity(x.rows(), x.cols());
    }
};

// Over a field they form a group, if we exclude singular matrices (a precondition).
template <SquareMatrix M>
    requires Field<typename M::Scalar>
struct inverse_operation<M, std::multiplies<>> {
    static M of(const M& x) { return x.inverse(); }
};

// Eigen declares operator* for matrices of any shape and checks the sizes with a
// static_assert inside it, so a 2×3 matrix looks like a semigroup to the concept.
template <EigenMatrix M>
    requires(!SquareMatrix<M>)
inline constexpr bool disable_semigroup<M, std::multiplies<>> = true;

// ---- Real and complex exponents --------------------------------------------

namespace detail {

// std::complex over a floating-point type: integer matrices are promoted to complex<double>.
template <class S>
using complex_of = std::complex<promote_t<typename Eigen::NumTraits<S>::Real>>;

}  // namespace detail

// Aᵗ = V · diag(λᵢᵗ) · V⁻¹ for a diagonalizable A, principal branch. Every λᵢᵗ is
// gp::pow(λᵢ, t) from the complex layer. The result is always a complex matrix,
// because a real matrix can have a complex power (the square root of a reflection).
template <SquareMatrix M, class E>
    requires Real<E> || Complex<E>
auto pow(const M& a, const E& t) {
    using C = detail::complex_of<typename M::Scalar>;
    using R = typename C::value_type;
    using CM = Eigen::Matrix<C, M::RowsAtCompileTime, M::ColsAtCompileTime>;
    using CV = Eigen::Matrix<C, M::RowsAtCompileTime, 1>;

    GP_EXPECTS(a.rows() == a.cols());
    const Eigen::ComplexEigenSolver<CM> es(a.template cast<C>());
    if (es.info() != Eigen::Success) throw std::domain_error("gp::pow: eigendecomposition failed");
    const Eigen::PartialPivLU<CM> lu(es.eigenvectors());
    if (lu.rcond() < std::sqrt(std::numeric_limits<R>::epsilon()))
        throw std::domain_error("gp::pow: matrix is not diagonalizable");

    CV lambda = es.eigenvalues();
    for (C& l : lambda) {
        // Rounding noise such as Im λ = −0.0 must not push a negative λ across the branch cut.
        if (std::abs(l.imag()) <= R(64) * std::numeric_limits<R>::epsilon() * std::abs(l)) l = C(l.real(), R(0));
        l = gp::pow(l, t);
    }
    return CM(es.eigenvectors() * lambda.asDiagonal() * lu.inverse());
}

// Expressions (A * B, A.transpose(), 2 * A) are evaluated first.
template <class D, class E>
    requires(!EigenMatrix<D>)
auto pow(const Eigen::MatrixBase<D>& x, const E& e) -> decltype(gp::pow(x.eval(), e)) {
    return gp::pow(x.eval(), e);
}

}  // namespace gp
