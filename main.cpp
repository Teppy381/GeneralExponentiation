// A tour of gp::power and gp::pow, one layer at a time. The [tag] on every row is
// the strongest algebraic structure the concepts found: it picks the power() overload.
#include <gp/gp.hpp>

#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <expected>
#include <format>
#include <functional>
#include <numbers>
#include <print>
#include <string>
#include <string_view>
#include <system_error>
#include <type_traits>

// ---------------- SHOWCASE OF CONSTEXPR FUNCTIONALITY ----------------

// -- Retroactive modeling: std::string with concatenation is a monoid, e = "" --
template <> struct gp::identity_element<std::string, std::plus<>> {
    static constexpr std::string of(const std::string &) { return {}; }
};

// -- A semigroup the library has never heard of: 2×2 matrices as std::array --
using mat2 = std::array<std::uint64_t, 4>; // [a b; c d]

constexpr auto mat2_mul = [](const mat2 &x, const mat2 &y) -> mat2 {
    return {x[0] * y[0] + x[1] * y[2], x[0] * y[1] + x[1] * y[3],
            x[2] * y[0] + x[3] * y[2], x[2] * y[1] + x[3] * y[3]};
};
using mat2_mul_t = std::remove_const_t<decltype(mat2_mul)>;

// [1 1; 1 0]ⁿ = [F(n+1) F(n); F(n) F(n−1)]. A lambda knows no identity element, so n > 0.
constexpr std::uint64_t fib(unsigned n) {
    return gp::power(mat2{1, 1, 1, 0}, n, mat2_mul)[1];
}

// Whatever needs no <cmath> is checked at compile time.
static_assert(fib(90) == 2880067194370816120ULL);
static_assert(gp::power(7, 13, std::plus<>{}) == 91);
static_assert(gp::power(std::string("ab"), 3, std::plus<>{}) == "ababab");
static_assert(gp::pow(3, 13) == 1594323);
static_assert(gp::pow(2.0, -3) == 0.125);
static_assert(!gp::checked_pow(2, 31));

namespace {

using gp::structure_of;
using std::numbers::pi;

using cd = std::complex<double>;
using Matrix2l = Eigen::Matrix<std::int64_t, 2, 2>;

const cd i{0, 1};

Eigen::Matrix2d rot(double a) {
    return Eigen::Matrix2d{{std::cos(a), -std::sin(a)}, {std::sin(a), std::cos(a)}};
}

std::string show(std::integral auto x) { return std::format("{}", x); }
// Prints −0 as 0.
std::string show(double x) { return std::format("{}", x == 0 ? 0.0 : x); }
std::string show(cd z) {
    return std::format("{}{:+}i", z.real() == 0 ? 0.0 : z.real(),
                       z.imag() == 0 ? 0.0 : z.imag());
}

template <class T> std::string show(const std::expected<T, std::errc> &r) {
    return r ? show(*r)
             : std::format("error: {}", std::make_error_code(r.error()).message());
}

// Rounding noise → 0.
double chop(double x) { return std::abs(x) < 1e-12 ? 0.0 : x; }

std::string entry(std::integral auto x) { return std::format("{}", x); }
std::string entry(double x) { return std::format("{:.6g}", chop(x)); }
std::string entry(cd z) {
    const double re = chop(z.real()), im = chop(z.imag());
    if (im == 0)
        return entry(re);
    if (re == 0)
        return std::format("{:.6g}i", im);
    return std::format("{:.6g}{:+.6g}i", re, im);
}

template <class D> std::string show(const Eigen::MatrixBase<D> &expr) {
    const typename D::PlainObject m = expr;
    std::string s;
    for (Eigen::Index r = 0; r < m.rows(); ++r)
        for (Eigen::Index c = 0; c < m.cols(); ++c)
            s += std::format("{}{}", c ? ", " : r ? "; " : "[", entry(m(r, c)));
    return s + "]";
}

void section(std::string_view title) { std::println("\n{}", title); }

void row(std::string_view expr, std::string_view value, std::string_view structure = "",
         std::string_view note = "") {
    std::string tail = structure.empty() ? std::string{} : std::format("[{}]", structure);
    if (!note.empty())
        tail += std::format("{}{}", tail.empty() ? "" : " ", note);
    std::string line = std::format("  {:<31} = {:<26} {}", expr, value, tail);
    line.erase(line.find_last_not_of(' ') + 1);
    std::println("{}", line);
}

void one_algorithm() {
    section("0. One algorithm, many structures: power(x, n, op)");

    row("power(7, 13, plus)", show(gp::power(7, 13, std::plus<>{})),
        structure_of<int, std::plus<>>(), "Egyptian multiplication: 13·7");
    row("power(7, -3, plus)", show(gp::power(7, -3, std::plus<>{})),
        structure_of<int, std::plus<>>(), "n < 0 through the inverse −x");
    row("power(\"ab\", 3, plus)", gp::power(std::string("ab"), 3, std::plus<>{}),
        structure_of<std::string, std::plus<>>(),
        "not commutative; e = \"\" is declared in main.cpp");
    row("power(Q, 90, lambda)[1]", show(fib(90)), structure_of<mat2, mat2_mul_t>(),
        "fib(90) on std::array: no e, so n > 0");
}

void integers() {
    section("1. Integers: a monoid under ×, so n >= 0");

    row("pow(3, 13)", show(gp::pow(3, 13)), structure_of<int>());
    row("pow(2LL, 62)", show(gp::pow(2LL, 62)), structure_of<long long>());
    const auto checked = structure_of<std::expected<long long, std::errc>,
                                      gp::checked_multiplies<long long>>();
    row("checked_pow(10LL, 18)", show(gp::checked_pow(10LL, 18)), checked,
        "same power(), std::expected monoid");
    row("checked_pow(10LL, 19)", show(gp::checked_pow(10LL, 19)), checked,
        "overflow is a value, not UB");
}

void reals() {
    section("2. Reals: a group under × (0 aside); real exponents build on integer ones");

    row("pow(2.0, -3)", show(gp::pow(2.0, -3)), structure_of<double>());
    row("pow(2.0, 0.5)", show(gp::pow(2.0, 0.5)), "",
        std::format("std::sqrt(2) = {}", std::sqrt(2.0)));
    row("pow(-8.0, 3.0)", show(gp::pow(-8.0, 3.0)), "", "integral y → the exact path");
    row("pow(-8.0, 1.0 / 3)", show(gp::pow(-8.0, 1.0 / 3)), "", "not a real number");
    row("pow(10, 0.5)", show(gp::pow(10, 0.5)), "", "an int base is promoted to double");
}

void complex_numbers() {
    section("3. Complex numbers: build on the reals, |z|ᵃ comes from layer 2");

    row("pow(i, 2)", show(gp::pow(i, 2)), structure_of<cd>());
    row("pow(1 + i, 8.0)", show(gp::pow(1.0 + i, 8.0)), "", "integral y → exact power()");
    row("std::pow(1 + i, 8.0)", show(std::pow(1.0 + i, 8.0)), "",
        "polar form, for comparison");
    row("pow(e, iπ)", show(gp::pow(std::numbers::e, i * pi)), "", "Euler: ≈ −1");
    row("pow(i, i)", show(gp::pow(i, i)), "", "= e^(−π/2), a real number");
}

void matrices() {
    section("4. Matrices: the same power(), Identity() as e");

    const Matrix2l Q{{1, 1}, {1, 0}};
    row("pow(Q, 90), Q = [1 1; 1 0]", show(gp::pow(Q, 90)), structure_of<Matrix2l>(),
        "Fibonacci in int64");
    row("pow(R(30°), 3)", show(gp::pow(rot(pi / 6), 3)), structure_of<Eigen::Matrix2d>(),
        "= R(90°)");
    const Eigen::Matrix2d A{{4, 7}, {2, 6}};
    row("pow(A, -3) · pow(A, 3)", show(gp::pow(A, -3) * gp::pow(A, 3)),
        structure_of<Eigen::Matrix2d>(), "n < 0 through A⁻¹");
    Eigen::MatrixXd M(3, 3);
    M << 2, 0, 1, 1, 3, 0, 0, 1, 4;
    row("pow(MatrixXd(3, 3), 0)", show(gp::pow(M, 0)), structure_of<Eigen::MatrixXd>(),
        "e sized by the sample");
    row("Matrix2i | Matrix2d | Matrix2cd",
        std::format("{} | {} | {}", structure_of<Eigen::Matrix2i>(),
                    structure_of<Eigen::Matrix2d>(), structure_of<Eigen::Matrix2cd>()),
        "", "A⁻¹ exists only over a field");
}

void complex_matrices() {
    section("5. Complex matrices: matrices over the field ℂ; Aᵗ builds on ℂ^t");

    const Eigen::Matrix2cd sx{{0, 1}, {1, 0}}, sy{{0, -i}, {i, 0}};
    row("pow(σy, 2)", show(gp::pow(sy, 2)), structure_of<Eigen::Matrix2cd>(),
        "Pauli matrices square to I");
    const double theta = pi / 3;
    const long n = 1L << 20;
    const Eigen::Matrix2cd step =
        Eigen::Matrix2cd::Identity() + (i * theta / static_cast<double>(n)) * sx;
    row("pow(I + iθσx/n, n), n = 2²⁰", show(gp::pow(step, n)),
        structure_of<Eigen::Matrix2cd>(), "20 squarings, θ = π/3");
    row("cos θ·I + i·sin θ·σx",
        show(std::cos(theta) * Eigen::Matrix2cd::Identity() + i * std::sin(theta) * sx),
        "", "= e^(iθσx), the limit");
    const Eigen::Matrix2d X{{0, 1}, {1, 0}};
    const auto sqrt_not = gp::pow(X, 0.5);
    row("pow(σx, 0.5)", show(sqrt_not), "", "√NOT: a real matrix with a complex root");
    row("pow(pow(σx, 0.5), 2)", show(gp::pow(sqrt_not, 2)),
        structure_of<Eigen::Matrix2cd>());
    row("pow(R(90°), 0.5)", show(gp::pow(rot(pi / 2), 0.5)), "", "= R(45°)");
}

} // namespace

int main() {
    one_algorithm();
    integers();
    reals();
    complex_numbers();
    matrices();
    complex_matrices();

    // Rejected at compile time:
    //   gp::pow(Eigen::Matrix<double, 2, 3>{}, 2);  // 2×3 is not a semigroup
    //   gp::pow(std::string("ab"), 2);              // no operator* for std::string
    //   static_assert(fib(0) == 0);                 // n > 0: a lambda knows no identity
    //   static_assert(gp::pow(2, -1) == 0);         // n >= 0: ℤ is a monoid, not a group
}
