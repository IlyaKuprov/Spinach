# kernel/utilities/clebsch_gordan.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/clebsch_gordan.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/clebsch_gordan.m)

## Purpose

Computes the Clebsch-Gordan coefficient: the coefficient in front of the `Y(L,M)` spherical harmonic in the expansion of the product of `Y(L1,M1)` and `Y(L2,M2)` spherical harmonics. In the more general sense, the coefficient refers to the expansion coefficient of the `|L,M>` angular momentum or spin state in the product basis of `|L1,M1>|L2,M2>` states.

## Behaviour

- Syntax: `cg=clebsch_gordan(L,M,L1,M1,L2,M2)`.
- The notation is matched to Varshalovich, Section 8.2.1, with the mapping `c=L`, `gam=M`, `a=L1`, `alp=M1`, `b=L2`, `bet=M2`.
- The function runs three stages of zero tests before computing anything:
  - Stage I checks that `a+alp>=0`, `a-alp>=0`, `b+bet>=0`, `b-bet>=0`, `c+gam>=0`, `c-gam>=0`, and `gam==alp+bet`.
  - Stage II checks the triangle conditions `a+b-c>=0`, `a-b+c>=0`, `-a+b+c>=0`, and `a+b+c+1>=0`.
  - Stage III checks that the summation range is non-empty, with lower limit `max([alp-a, b+gam-a, 0])` and upper limit `min([c+b+alp, c+b-a, c+gam])`.
- If any test fails, the function returns zero for inadmissible index combinations.
- If all tests pass, the coefficient is computed using arbitrary-precision Java arithmetic:
  - A look-up table of factorials is built with `java.math.BigInteger`, with `gamma(k)` holding `(k-1)!` for `k` from 3 to `a+b+c+2`.
  - The prefactor numerator is `(2c+1) * gamma(a-b+c+1) * gamma(-a+b+c+1) * gamma(c+gam+1) * gamma(c-gam+1) * gamma(a+b-c+1)`.
  - The prefactor denominator is `gamma(a+b+c+2) * gamma(a+alp+1) * gamma(a-alp+1) * gamma(b+bet+1) * gamma(b-bet+1)`.
  - The summation runs over `z` from the lower to the upper limit, with each term being `(-1)^(b+bet+z) * gamma(c+b+alp-z+1) * gamma(a-alp+z+1)` divided by `gamma(z+1) * gamma(c-a+b-z+1) * gamma(c+gam-z+1) * gamma(a-b-gam+z+1)`.
  - Divisions use `java.math.BigDecimal` with `java.math.RoundingMode.HALF_UP`, at a precision of `length(denom.toString)-length(numer.toString)+64` digits.
  - The final result is returned as `sign(z_sum) * sqrt(z_sum^2 * prefactor)` converted to double precision.
- Per the source header, the calculation in double-precision arithmetic is not trivial for high ranks; this function produces machine-precision answers up to about `L=1e4`. A faster implementation for low ranks is available in `cg_fast.m`.

## Inputs and outputs

**Inputs**

- `L, M, L1, M1, L2, M2` — integer or half-integer indices of the angular momentum or spin states. Each must be numeric, real, and a scalar, and each must satisfy `mod(2*x+1,1)==0` (i.e., be integer or half-integer). `L`, `L1`, and `L2` must not exceed `1e4`.

**Outputs**

- `cg` — floating-point (double precision) Clebsch-Gordan coefficient. Zero is returned for inadmissible index combinations.

**Errors**

- `'all arguments must be numeric.'` if any argument is not numeric.
- `'all arguments must be real.'` if any argument is not real.
- `'all arguments must have one element.'` if any argument is not a scalar.
- `'all arguments must be integer or half-integer.'` if any argument fails the integrality test.
- `'you must be joking.'` if `L`, `L1`, or `L2` exceeds `1e4`.

## References

- Varshalovich, D. A., Section 8.2.1 (notation reference cited in the source).
- Spinach Wiki: [clebsch_gordan.m](https://spindynamics.org/wiki/index.php?title=clebsch_gordan.m)
