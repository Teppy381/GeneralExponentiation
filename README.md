# GeneralExponentiation

Generic exponentiation in C++23. One algorithm, `gp::power(x, n, op)`, raises integers,
reals, complex numbers, matrices and complex matrices, because all it needs is an
associative operation. Concepts pick the overload, and each numeric layer is built on the
previous one. Header-only; needs Clang 20 or GCC 14, CMake 3.28 and Eigen 5, or newer.

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
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
./build/showcase
```

The build files set no compiler-specific flags. Warnings and sanitizers go on the command line.

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_COMPILER=clang++ \
      -DCMAKE_CXX_FLAGS="-Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion -Wshadow -Werror -fsanitize=address,undefined"
```

Eigen comes from the system package (`eigen3-devel` on Fedora). Without it, CMake fetches
the Eigen 5.0.1 headers. The tests need doctest and are skipped without it. Clang 18 is
too old: it cannot compile `std::expected` from libstdc++. GCC 13 compiles the headers, but
the showcase needs `<print>` from libstdc++ 14.

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

## Showcase

`./build/showcase` prints one row per call. The `[tag]` is the strongest
structure the concepts found for the type, and so the `power` overload that ran.

```
0. One algorithm, many structures: power(x, n, op)
  power(7, 13, plus)              = 91                         [Group] Egyptian multiplication: 13·7
  power(7, -3, plus)              = -21                        [Group] n < 0 through the inverse −x
  power("ab", 3, plus)            = ababab                     [Monoid] not commutative; e = "" is declared in main.cpp
  power(Q, 90, lambda)[1]         = 2880067194370816120        [Semigroup] fib(90) on std::array: no e, so n > 0

1. Integers: a monoid under ×, so n >= 0
  pow(3, 13)                      = 1594323                    [Monoid]
  pow(2LL, 62)                    = 4611686018427387904        [Monoid]
  checked_pow(10LL, 18)           = 1000000000000000000        [Monoid] same power(), std::expected monoid
  checked_pow(10LL, 19)           = error: Numerical result out of range [Monoid] overflow is a value, not UB

2. Reals: a group under × (0 aside); real exponents build on integer ones
  pow(2.0, -3)                    = 0.125                      [Group]
  pow(2.0, 0.5)                   = 1.414213562373095          std::sqrt(2) = 1.4142135623730951
  pow(-8.0, 3.0)                  = -512                       integral y → the exact path
  pow(-8.0, 1.0 / 3)              = nan                        not a real number
  pow(10, 0.5)                    = 3.1622776601683795         an int base is promoted to double

3. Complex numbers: build on the reals, |z|ᵃ comes from layer 2
  pow(i, 2)                       = -1+0i                      [Group]
  pow(1 + i, 8.0)                 = 16+0i                      integral y → exact power()
  std::pow(1 + i, 8.0)            = 15.999999999999998-3.9188697572715295e-15i polar form, for comparison
  pow(e, iπ)                      = -1+1.2246467991473532e-16i Euler: ≈ −1
  pow(i, i)                       = 0.20787957635076193+0i     = e^(−π/2), a real number

4. Matrices: the same power(), Identity() as e
  pow(Q, 90), Q = [1 1; 1 0]      = [4660046610375530309, 2880067194370816120; 2880067194370816120, 1779979416004714189] [Monoid] Fibonacci in int64
  pow(R(30°), 3)                  = [0, -1; 1, 0]              [Group] = R(90°)
  pow(A, -3) · pow(A, 3)          = [1, 0; 0, 1]               [Group] n < 0 through A⁻¹
  pow(MatrixXd(3, 3), 0)          = [1, 0, 0; 0, 1, 0; 0, 0, 1] [Group] e sized by the sample
  Matrix2i | Matrix2d | Matrix2cd = Monoid | Group | Group     A⁻¹ exists only over a field

5. Complex matrices: matrices over the field ℂ; Aᵗ builds on ℂ^t
  pow(σy, 2)                      = [1, 0; 0, 1]               [Group] Pauli matrices square to I
  pow(I + iθσx/n, n), n = 2²⁰     = [0.5, 0.866026i; 0.866026i, 0.5] [Group] 20 squarings, θ = π/3
  cos θ·I + i·sin θ·σx            = [0.5, 0.866025i; 0.866025i, 0.5] = e^(iθσx), the limit
  pow(σx, 0.5)                    = [0.5+0.5i, 0.5-0.5i; 0.5-0.5i, 0.5+0.5i] √NOT: a real matrix with a complex root
  pow(pow(σx, 0.5), 2)            = [0, 1; 1, 0]               [Group]
  pow(R(90°), 0.5)                = [0.707107, -0.707107; 0.707107, 0.707107] = R(45°)
```

## Caveats

- Include `gp/gp.hpp` only. A trait specialization must be visible wherever a concept is
  checked; otherwise satisfaction differs between translation units, which makes the
  program ill-formed with no diagnostic required.
- `ℝ^ℤ` by squaring has a worst-case relative error of about `(n − 1)·ε`, like a naive
  product, while libm's `pow` stays within an ulp. That is the price of one algorithm for
  numbers and matrices alike.
- `ℝ^ℝ` inherits that error through `xⁿ`, and `ℂ^ℝ`, `ℂ^ℂ` and `Aᵗ` inherit it through
  `|z|ᵃ`. Against `exp(y·ln x)`, whose error grows as `|y·ln x|·ε`, the split is a few
  times more accurate when `x` is far from 1 and much worse when it is close: for
  `x ≈ 1 + 10⁻⁷` and `y ≈ 10⁷` only about 9 digits are correct. It is also up to 10 times
  slower than `std::pow`. The layer shows how real exponents reduce to integer ones; it
  does not compete with libm.
- `Aᵗ` requires a diagonalizable `A`: the conditioning of the eigenbasis is checked, and
  `std::domain_error` is thrown otherwise. The result is always a complex matrix.
- With Clang and libstdc++, `std::complex` arithmetic cannot be used in constant
  expressions, so complex results are checked at run time.
- A negative power of an `int` or of a singular matrix is a precondition violation: a
  compile error in a constant expression, an abort at run time.
