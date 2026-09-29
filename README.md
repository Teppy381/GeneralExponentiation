# GeneralExponentiation

Generic exponentiation in C++23. One algorithm, `gp::power(x, n, op)`, raises integers,
reals, complex numbers, matrices and complex matrices, because all it needs is an
associative operation. Concepts pick the overload, and each numeric layer is built on the
previous one. Header-only; tested with Clang 22, libstdc++ 16, CMake 4.3 and Eigen 5.

```cpp
#include <gp/gp.hpp>

gp::pow(3, 13);                      // 1594323, in constant expressions too
gp::pow(2.0, -3);                    // 0.125
gp::pow(2.0, 0.5);                   // 1.414213562373095
gp::pow(std::complex{0.0, 1.0}, 2);  // -1+0i
gp::pow(a, -3);                      // a⁻³ for an Eigen::Matrix2d a
gp::pow(pauli_x, 0.5);               // √NOT, a complex matrix
gp::power(7, 13, std::plus<>{});     // 91: with + the same algorithm multiplies
gp::checked_pow(10LL, 19);           // std::unexpected(std::errc::result_out_of_range)
```

## The idea

`xⁿ` needs nothing but associativity: `x⁸ = ((x²)²)²` is `x·x·…·x` with the parentheses
moved. Every extra axiom widens the set of exponents that make sense:

| Concept | Adds | Examples | `xⁿ` for |
|---|---|---|---|
| `Semigroup<T, Op>` | associative `op` | a lambda multiplying `std::array`s | `n > 0` |
| `Monoid<T, Op>` | identity `e` | `(int, ×, 1)`, `(std::string, +, "")`, `(Matrix2i, ×, I)` | `n ≥ 0` |
| `Group<T, Op>` | inverse `x⁻¹` | `(double, ×)`, `(int, +)`, `(Matrix2d, ×)` | any `n` |

Each concept refines the previous one, so the three `power` overloads are ordered by
subsumption: the compiler picks the strongest one that applies. The algorithm is
Stepanov's (*From Mathematics to Generic Programming*, ch. 7): `⌊log₂ n⌋` squarings plus
`popcount(n) − 1` multiplications, and no squaring past the top bit. So no intermediate
value exceeds the result, and that is what makes `checked_pow` exact.

Types are adapted without being touched, by specializing traits:

```cpp
template <>
struct gp::identity_element<std::string, std::plus<>> {
    static constexpr std::string of(const std::string&) { return {}; }
};

static_assert(gp::power(std::string("ab"), 3, std::plus<>{}) == "ababab");
```

`of` takes a sample because an identity can depend on a run-time size (`MatrixXd`).
`gp::inverse_operation<T, Op>` works the same way. `gp::disable_semigroup<T, Op>` is the
opt-out for types that multiply syntactically but not mathematically: Eigen declares
`operator*` for any two matrices and rejects a 2×3 · 2×3 product only inside its body.

## Layers

```
power(x, n, op)                         Stepanov, O(log n)
 ├─ ℤ^ℤ, ℝ^ℤ, ℂ^ℤ, Mat^ℤ, Mat<ℂ>^ℤ      the same power(); only the traits differ
 └─ ℝ^ℝ = xⁿ · e^{f·ln x}                on ℝ^ℤ, with n = trunc y, f = y − n
     └─ ℂ^ℝ, ℂ^ℂ = |z|ᵃ · e^{−bθ} · cis(b·ln|z| + aθ)    on ℝ^ℝ
         └─ Aᵗ = V · diag(λᵢᵗ) · V⁻¹      on ℂ^t, for real and complex matrices
```

| Call | Overload | Algorithm |
|---|---|---|
| `pow(3, 13)`, `pow(2.0, -3)`, `pow(z, 8)`, `pow(a, 90)` | `pow(const T&, Integer)` | `power(x, n, *)` |
| `checked_pow(10LL, 19)` | `checked_pow(Integer, Integer)` | `power` over `std::expected` |
| `pow(2.0, 0.5)`, `pow(2, 0.5)` | `pow(Arithmetic, Real)` | an integral `y` takes the exact path |
| `pow(z, 0.5)`, `pow(z, w)`, `pow(e, iπ)` | `pow(Complex, Real)`, `pow(Complex, Complex)`, `pow(Arithmetic, Complex)` | principal branch |
| `pow(a, 0.5)` | `pow(SquareMatrix, Real \| Complex)` | `Eigen::ComplexEigenSolver` |
| `pow(a * b, 3)` | `pow(const Eigen::MatrixBase<D>&, E)` | `.eval()`, then one of the above |

## Build

```sh
cmake --preset clang-debug         # Debug with ASan + UBSan and -Werror
cmake --build --preset clang-debug
ctest --preset clang-debug
./build/clang-debug/showcase
```

Eigen comes from the system package (`eigen3-devel` on Fedora). Without it, CMake fetches
the Eigen 5.0.1 headers. The tests need doctest and are skipped without it.

```
include/gp/
  concepts.hpp   domain and algebraic concepts, customization points
  power.hpp      power(x, n, op) ×3, pow(x, n)
  scalar.hpp     ℤ and ℝ: traits, checked_pow, real exponents
  complex.hpp    ℂ: traits, complex exponents
  matrix.hpp     the Eigen adapter, Aᵗ
  gp.hpp         the umbrella header
main.cpp         a tour of the layers
tests/           static_asserts and doctest cases
```

## Caveats

- Include `gp/gp.hpp` only. A trait specialization must be visible wherever a concept is
  checked; otherwise satisfaction differs between translation units, which makes the
  program ill-formed with no diagnostic required.
- `ℝ^ℤ` by squaring has a worst-case relative error of about `(n − 1)·ε`, like a naive
  product, while libm's `pow` stays within an ulp. That is the price of one algorithm for
  numbers and matrices alike.
- `Aᵗ` requires a diagonalizable `A`: the conditioning of the eigenbasis is checked, and
  `std::domain_error` is thrown otherwise. The result is always a complex matrix.
- With Clang and libstdc++, `std::complex` arithmetic cannot be used in constant
  expressions, so complex results are checked at run time.
- A negative power of an `int` or of a singular matrix is a precondition violation: a
  compile error in a constant expression, an abort at run time.
