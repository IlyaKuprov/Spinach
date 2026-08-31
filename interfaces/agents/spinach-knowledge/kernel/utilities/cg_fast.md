# kernel/utilities/cg_fast.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/cg_fast.m`
- Signature: `cg=cg_fast(L,M,L1,M1,L2,M2)`
- Total lines: 111

## Purpose

Clebsch-Gordan coefficient: the coefficient in front of Y(L,M) spheri- cal harmonic in the expansion of the product of Y(L1,M1) and Y(L2,M2) spherical harmonics. In the more general sense, the coefficient refers to the expansion coefficient of |L,M> angular momentum or spin state in the product basis of |L1,M1>|L2,M2> states. Syntax: cg=cg_fast(L,M,L1,M1,L2,M2)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 35-36: Check consistency; implemented by `grumble(L,M,L1,M1,L2,M2)`.
- Lines 38-39: Match the notation to Varshalovich, Section 8.2.1; implemented by `cg=0; c=L; gam=M; a=L1; alp=M1; b=L2; bet=M2`.
- Lines 41-45: Run zero tests (Stage I); implemented by `prefactor_is_nonzero=(a+alp>=0)&&(a-alp>=0)&& (b+bet>=0)&&(b-bet>=0)&& (c+gam>=0)&&(c-gam>=0)&& (gam==alp+bet)`.
- Lines 47-48: Proceed if appropriate; implemented by `if prefactor_is_nonzero`.
- Lines 50-52: Run zero tests (Stage II); implemented by `delta_is_nonzero=(a+b-c>=0)&&(a-b+c>=0)&& (-a+b+c>=0)&&(a+b+c+1>=0)`.
- Lines 54-55: Proceed if appropriate; implemented by `if delta_is_nonzero`.
- Lines 57-58: Run zero tests (Stage III); implemented by `lower_sum_limit=max([alp-a, b+gam-a, 0])`.
- Lines 62-63: Proceed if appropriate; implemented by `if sum_is_nonzero`.
- Lines 65-66: Compute Equation 8.2.1(5) using log-factorials; implemented by `for z=lower_sum_limit:upper_sum_limit`.

### Control flow inferred from the code

- Line 48: conditional branch on `prefactor_is_nonzero`.
- Line 55: conditional branch on `delta_is_nonzero`.
- Line 63: conditional branch on `sum_is_nonzero`.
- Line 66: `for` loop over `z=lower_sum_limit:upper_sum_limit`.

### Key state/data transformations

- Lines 39: computes `cg` using `cg=0; c=L; gam=M; a=L1; alp=M1; b=L2; bet=M2`.
- Lines 58: computes `lower_sum_limit` using `lower_sum_limit=max([alp-a, b+gam-a, 0])`.
- Lines 59: computes `upper_sum_limit` using `upper_sum_limit=min([c+b+alp, c+b-a, c+gam])`.

### Local helper functions

- Line 85: `grumble()` — `function grumble(L,M,L1,M1,L2,M2)`.
  - Representative operation: `if (~isnumeric(L))||(~isnumeric(L1))||(~isnumeric(L2))|| (~isnumeric(M))||(~isnumeric(M1))||(~isnumeric(M2))`.
  - Representative operation: `(~isnumeric(M))||(~isnumeric(M1))||(~isnumeric(M2))`.

## Parameters / inputs

- L,M,L1,M1,L2,M2 -integer or half-integer indices of
- the angular momentum or spin states

## Outputs

- cg -floating-point (double precision)
- Clebsch-Gordan coefficient
- Note: only some combinations of L,M,L1,M1,L2,M2 are allowed by the pro-
- perties of spherical harmonics and spin states. If inadmissible
- indices are supplied, zero is returned.
- Note: CG coefficient calculation in double-precision arithmetic is not
- a trivial matter for high ranks. This function produces fast ans-
- wers with an accuracy of about 1e-3 up to about L=20. A slower
- machine precision implementation for higher ranks is available
- in clebsch_gordan.m function.

## Implementation structure

- Clebsch-Gordan coefficient: the coefficient in front of Y(L,M) spheri-
- cal harmonic in the expansion of the product of Y(L1,M1) and Y(L2,M2)
- spherical harmonics. In the more general sense, the coefficient refers
- to the expansion coefficient of |L,M> angular momentum or spin state
- in the product basis of |L1,M1>|L2,M2> states. Syntax:
- cg=cg_fast(L,M,L1,M1,L2,M2)
- L,M,L1,M1,L2,M2 -integer or half-integer indices of
- the angular momentum or spin states
- cg -floating-point (double precision)
- Clebsch-Gordan coefficient
- Note: only some combinations of L,M,L1,M1,L2,M2 are allowed by the pro-
- perties of spherical harmonics and spin states. If inadmissible

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `logfactorial()`, `clebsch_gordan()`.
