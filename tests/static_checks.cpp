// Compile-time tests: this file compiles only if they pass.
#include <gp/gp.hpp>

#include <complex>
#include <functional>
#include <limits>
#include <string>

namespace {

using cd = std::complex<double>;

template <class T, class E>
concept can_pow = requires(const T &x, const E &e) { gp::pow(x, e); };

// ---- The structures the concepts see ----
static_assert(gp::Monoid<int, std::multiplies<>> && !gp::Group<int, std::multiplies<>>);
static_assert(gp::Group<double, std::multiplies<>> && gp::Group<cd, std::multiplies<>>);
static_assert(gp::Group<unsigned, std::plus<>>);
static_assert(gp::Monoid<Eigen::Matrix2i, std::multiplies<>> &&
              !gp::Group<Eigen::Matrix2i, std::multiplies<>>);
static_assert(gp::Group<Eigen::Matrix2d, std::multiplies<>> &&
              gp::Group<Eigen::MatrixXcd, std::multiplies<>>);
static_assert(!gp::Semigroup<Eigen::Matrix<double, 2, 3>, std::multiplies<>>);
static_assert(!gp::Semigroup<std::string, std::multiplies<>>);

// ---- The overload set ----
static_assert(can_pow<int, int> && can_pow<double, double> && can_pow<int, float>);
static_assert(can_pow<cd, int> && can_pow<cd, double> && can_pow<cd, cd> &&
              can_pow<double, cd>);
static_assert(can_pow<Eigen::MatrixXd, int> && can_pow<Eigen::Matrix2i, double> &&
              can_pow<Eigen::Matrix2cd, cd>);
static_assert(!can_pow<Eigen::Matrix<double, 2, 3>, int> &&
              !can_pow<Eigen::Matrix<double, 2, 3>, double>);
static_assert(!can_pow<int, bool> && !can_pow<std::string, int>);

// ---- Values ----
static_assert(gp::pow(3, 13) == 1594323);
static_assert(gp::pow(-2, 0) == 1);
static_assert(gp::pow(2.0, -3) == 0.125);
static_assert(gp::power(7, -3, std::plus<>{}) == -21);
static_assert(gp::checked_pow(-2, 31) == std::numeric_limits<int>::min()); // fits exactly
static_assert(!gp::checked_pow(2, 31) && !gp::checked_pow(3LL, 40));

// ---- checked_pow at the edges of the type ----
constexpr int int_min = std::numeric_limits<int>::min();
static_assert(gp::checked_pow(int_min, 1) == int_min && !gp::checked_pow(int_min, 2));
static_assert(!gp::checked_multiplies<int>{}(int_min, -1)); // |min| has no positive twin
static_assert(gp::checked_multiplies<int>{}(int_min, 1) == int_min);
static_assert(gp::checked_pow(-2, 32).error() == std::errc::result_out_of_range);
static_assert(gp::checked_pow(0, 5) == 0 && gp::checked_pow(-1, 5) == -1);
static_assert(gp::checked_pow(2U, 31) == 2147483648U && !gp::checked_pow(2U, 32));
static_assert(gp::checked_pow(short(-2), 15) == -32768 && !gp::checked_pow(short(2), 15));
static_assert(gp::checked_pow(10LL, 18) == 1'000'000'000'000'000'000 &&
              !gp::checked_pow(10LL, 19));

// ---- O(log n): ⌊log₂ n⌋ squarings + popcount(n) − 1 multiplications ----
struct counting_plus {
    int *count;
    constexpr int operator()(int a, int b) const { return ++*count, a + b; }
};

constexpr int operations(int n) {
    int count = 0;
    gp::power(1, n, counting_plus{&count});
    return count;
}

static_assert(operations(1'000'000) == 19 + 6);
static_assert(operations(1 << 20) == 20);

} // namespace
