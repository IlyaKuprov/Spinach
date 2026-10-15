# kernel/utilities/clebsch_gordan.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/clebsch_gordan.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/clebsch_gordan.m)

## Purpose

Computes the Clebsch-Gordan coefficient: the coefficient in front of the `Y(L,M)` spherical harmonic in the expansion of the product of `Y(L1,M1)` and `Y(L2,M2)` spherical harmonics. In the more general sense, the coefficient refers to the expansion coefficient of the `|L,M>` angular momentum or spin state in the product basis of `|L1,M1>|L2,M2>` states.

## Behaviour and numerical method

The three preliminary zero-test stages are retained: projection bounds and
`M=M1+M2`, triangle inequalities, and the legacy nonempty summation-range test.
A failed stage returns zero. After these retained stages, incompatible total-rank
or rank/projection integer–half-integer parities also return zero before the
recurrence is constructed. The notation follows Varshalovich, Section 8.2.1.

Admissible coefficients are evaluated using the Schulten–Gordon magnetic-index
three-term recurrence with Luscombe–Luban ratio propagation in the monotone tails.
Cyclic column permutations select the shortest family without an additional
permutation sign. Magnetic recurrence avoids the singular zero-angular-momentum
endpoint of the alternative angular-momentum recurrence. Neighbouring interior
values are fitted, the full family is normalised by 3j orthogonality, and its
endpoint phase is combined with the Condon–Shortley CG conversion.

Coefficients, ratios, and recurrence amplitudes retain roundoff corrections in a
second native double; this is fixed compensated arithmetic, not an arbitrary-
precision runtime. The final coefficient is an ordinary double. There is no Symbolic Math Toolbox dependency. No factorial or alternating factorial sum
is evaluated. Small nonzero tails are reconstructed before the final conversion,
without squaring the target coefficient. FP64 output range and subnormal rounding
still apply; no universal correctly-rounded guarantee is made.

This is the shared coefficient implementation for all ranks. The two shipped Mathematica comparison
tables and `tests/kernel/test_wigner_angular_suite.m` cover numerical regressions;
tiny nonzero coefficients require relative as well as absolute comparisons.

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
- Schulten and Gordon, J. Math. Phys. 16, 1961 (1975), https://doi.org/10.1063/1.522426.
- Luscombe and Luban, Phys. Rev. E 57, 7274 (1998), https://doi.org/10.1103/PhysRevE.57.7274.
- Spinach Wiki: [clebsch_gordan.m](https://spindynamics.org/wiki/index.php?title=clebsch_gordan.m)
