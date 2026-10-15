# tests/kernel/test_wigner_angular_suite.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_wigner_angular_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_wigner_angular_suite.m)

## Purpose

Regression test suite for Spinach angular-momentum coefficient and spherical-function helpers. The suite checks Clebsch-Gordan coefficients, Wigner symbols, Wigner D matrices, and spherical harmonics against elementary exact and high-precision values, verifying that the angular-momentum helpers reproduce elementary exact quantum-mechanical coefficients.

## Behaviour

Each check compares a production calculation with an independent value or a
physical identity, recording explanatory failures without losing later checks.

Clebsch-Gordan checks:

- `clebsch_gordan(1,0,1/2,1/2,1/2,-1/2)` is compared to `1/sqrt(2)` (tolerances `1e-14`), because the `|1,0>` triplet contains `|alpha beta>` with coefficient `1/sqrt(2)`.
- `clebsch_gordan(0,0,1/2,1/2,1/2,-1/2)` is compared to `1/sqrt(2)` (tolerances `1e-14`), because the `|0,0>` singlet contains `|alpha beta>` with coefficient `1/sqrt(2)` in the Spinach phase convention.
- `clebsch_gordan(1,1,1/2,1/2,1/2,-1/2)` is compared to `0` (tolerances `1e-14`), confirming that forbidden projection combinations return zero.

Wigner symbol checks:

- `wigner_3j(1,0,1,0,0,0)` is compared to `-1/sqrt(3)` (tolerances `1e-14`), the elementary `(1 1 0; 0 0 0)` Wigner 3j symbol.
- `wigner_6j(0,0,0,0,0,0)` is compared to `1` (tolerances `1e-14`), the all-zero Wigner 6j symbol being unity.
- `wigner_6j(1,1,1,1,1,1)` is compared to `1/6` (tolerances `1e-14`), the Wigner 6j symbol with all angular momenta one.

Wigner D matrix checks:

- `wigner(1,0,pi/2,0)` is compared elementwise to the 3-by-3 reference matrix `[1/2 -1/sqrt(2) 1/2; 1/sqrt(2) 0 -1/sqrt(2); 1/2 1/sqrt(2) 1/2]` (tolerances `1e-14`), the Brink-Satchler closed form for the rank-one Wigner matrix at `beta=pi/2`.
- `wigner(2,0,0,0)` is compared to `eye(5)` (tolerances `1e-14`), because zero Euler angles give the identity Wigner D matrix.
- For `D=wigner(2,0.2,0.4,0.7)`, the product `D'*D` is compared to `eye(5)` (tolerances `1e-13`), verifying that Wigner D matrices are unitary rotation representations.

Spherical harmonics checks, evaluated on the angle grids `th=[0 pi/2 pi]` and `ph=[0 pi/3 pi/7]`:

- `spher_harmon(0,0,th,ph)` is compared to `ones(size(th))/sqrt(4*pi)` (tolerances `1e-14`), because `Y_0^0` is the constant `1/sqrt(4*pi)`.
- `spher_harmon(1,0,th,ph)` is compared to `sqrt(3/(4*pi))*cos(th)` (tolerances `1e-14`), because `Y_1^0` is `sqrt(3/(4*pi))*cos(theta)`.

## Inputs and outputs

Syntax:

```matlab
result=test_wigner_angular_suite()
```

The function takes no inputs. It returns `result`, a regression test result with explanatory messages accumulated across all `test_close` checks.

## References

1. Spinach source file: [tests/kernel/test_wigner_angular_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_wigner_angular_suite.m)

## Clebsch–Gordan recurrence coverage

The CG checks include exact-input 100-digit Mathematica reference values across
integer and half-integer magnetic recurrence regimes, representable tails near
`1e-300` compared relatively, a rank-10000 singlet, and an exact antisymmetry zero.
The two full 10000-row fundamentals comparison examples remain separate checks.
