# kernel/utilities/cg_fast.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/cg_fast.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/cg_fast.m)

## Purpose

Computes the Clebsch-Gordan coefficient `cg_fast(L,M,L1,M1,L2,M2)`, i.e. the coefficient in front of the `Y(L,M)` spherical harmonic in the expansion of the product of `Y(L1,M1)` and `Y(L2,M2)` spherical harmonics. In the more general sense, the coefficient is the expansion coefficient of the `|L,M>` angular momentum or spin state in the product basis of `|L1,M1>|L2,M2>` states.

## Behaviour

- Syntax: `cg=cg_fast(L,M,L1,M1,L2,M2)`.
- The notation is matched to Varshalovich, Section 8.2.1, and the coefficient is computed from Equation 8.2.1(5) using log-factorials.
- The function runs three stages of zero tests before summation:
  - Stage I: non-negativity of `a+alp`, `a-alp`, `b+bet`, `b-bet`, `c+gam`, `c-gam`, and the selection rule `gam == alp+bet`.
  - Stage II: triangle-inequality-type conditions `a+b-c >= 0`, `a-b+c >= 0`, `-a+b+c >= 0`, `a+b+c+1 >= 0`.
  - Stage III: the summation range check, with lower limit `max([alp-a, b+gam-a, 0])` and upper limit `min([c+b+alp, c+b-a, c+gam])`; the sum is only evaluated if the upper limit is not below the lower limit.
- If inadmissible index combinations are supplied, zero is returned.
- Accuracy is about `1e-3` up to about `L=20`; the function errors if any of `L`, `L1`, `L2` exceeds 20, directing the user to `clebsch_gordan()` instead, which is a slower machine-precision implementation for higher ranks.
- Input validation (via the internal `grumble` function) errors unless all six arguments are numeric, real, single-element, and integer or half-integer (checked via `mod(2*x+1,1)==0`).

## Inputs and outputs

**Inputs**

- `L, M, L1, M1, L2, M2` — integer or half-integer indices of the angular momentum or spin states.

**Outputs**

- `cg` — floating-point (double precision) Clebsch-Gordan coefficient.

## References

- D. A. Varshalovich, A. N. Moskalev, V. K. Khersonskii, *Quantum Theory of Angular Momentum*, Section 8.2.1, Equation 8.2.1(5).
- Spinach Wiki: [cg_fast.m](https://spindynamics.org/wiki/index.php?title=cg_fast.m)
