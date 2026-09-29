#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <gp/gp.hpp>

#include <unsupported/Eigen/MatrixFunctions>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <functional>
#include <limits>
#include <ranges>
#include <stdexcept>
#include <string>

namespace {

using cd = std::complex<double>;

// Equal within a relative 1e-14, with ±0, ±inf and NaN compared exactly.
bool same(double a, double b) {
    if (std::isnan(a) || std::isnan(b)) return std::isnan(a) && std::isnan(b);
    if (!std::isfinite(b) || b == 0) return a == b && std::signbit(a) == std::signbit(b);
    return std::abs(a - b) <= 1e-14 * std::abs(b);
}

}  // namespace

TEST_CASE("power agrees with the naive fold") {
    for (int n = 1; n <= 20; ++n) {
        CHECK(gp::power(3LL, n) == std::ranges::fold_left(std::views::repeat(3LL, n), 1LL, std::multiplies<>{}));
        const std::string s = "ab";
        CHECK(gp::power(s, n, std::plus<>{}) ==
              std::ranges::fold_left(std::views::repeat(s, n), std::string{}, std::plus<>{}));
    }
}

TEST_CASE("real exponents agree with std::pow, special values included") {
    const double inf = std::numeric_limits<double>::infinity(), nan = std::numeric_limits<double>::quiet_NaN();
    for (double x : {0.0, -0.0, 0.5, 1.0, -1.0, 2.0, 10.0, -2.0, inf, -inf, nan})
        for (double y : {0.0, 0.5, -0.5, 1.0, 3.0, -3.0, 2.5, -7.25, 1e20, -1e20, inf, -inf, nan}) {
            CAPTURE(x);
            CAPTURE(y);
            CHECK(same(gp::pow(x, y), std::pow(x, y)));
        }
}

TEST_CASE("complex exponents agree with std::pow") {
    for (cd z : {cd(1, 1), cd(-4, 0), cd(0, 1), cd(2, -3), cd(-1, -1)})
        for (cd w : {cd(2, 0), cd(0.5, 0), cd(-1.5, 0), cd(0, 1), cd(1, -2), cd(3, 0.5)}) {
            CAPTURE(z);
            CAPTURE(w);
            CHECK(std::abs(gp::pow(z, w) - std::pow(z, w)) <= 1e-12 * std::abs(std::pow(z, w)));
        }
    CHECK(gp::pow(cd(1, 1), 8.0) == cd(16, 0));  // an integral exponent is exact
}

TEST_CASE("matrix powers") {
    const Eigen::Matrix<std::int64_t, 2, 2> fib{{1, 1}, {1, 0}};
    CHECK(gp::pow(fib, 90)(0, 1) == 2880067194370816120);

    const Eigen::Matrix2d a{{4, 1}, {2, 3}};  // eigenvalues 2 and 5
    CHECK((gp::pow(a, -3) * gp::pow(a, 3)).isIdentity(1e-12));
    const Eigen::MatrixXd e = gp::pow(Eigen::MatrixXd::Random(3, 3), 0);  // an expression, sized at run time
    CHECK((e.rows() == 3 && e.isIdentity()));

    // Aᵗ agrees with Eigen's MatrixPower (Schur–Padé) and composes back.
    const Eigen::Matrix2cd root = gp::pow(a, 0.5);
    const Eigen::Matrix2d reference = a.pow(0.5);
    CHECK(root.isApprox(reference.cast<cd>(), 1e-12));
    CHECK(gp::pow(root, 2).isApprox(a.cast<cd>(), 1e-12));

    // √NOT = σx^½ = ½[1+i 1−i; 1−i 1+i] on the principal branch.
    const Eigen::Matrix2cd sqrt_not = gp::pow(Eigen::Matrix2d{{0, 1}, {1, 0}}, 0.5);
    const Eigen::Matrix2cd expected{{cd(0.5, 0.5), cd(0.5, -0.5)}, {cd(0.5, -0.5), cd(0.5, 0.5)}};
    CHECK(sqrt_not.isApprox(expected, 1e-12));

    // A Jordan block has no eigenbasis.
    CHECK_THROWS_AS(gp::pow(Eigen::Matrix2d{{1, 1}, {0, 1}}, 0.5), std::domain_error);
}
