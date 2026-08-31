# kernel/utilities/clebsch_gordan.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/clebsch_gordan.m`
- Signature: `cg=clebsch_gordan(L,M,L1,M1,L2,M2)`
- Total lines: 157

## Purpose

Clebsch-Gordan coefficient: the coefficient in front of Y(L,M) spheri- cal harmonic in the expansion of the product of Y(L1,M1) and Y(L2,M2) spherical harmonics. In the more general sense, the coefficient refers to the expansion coefficient of |L,M> angular momentum or spin state in the product basis of |L1,M1>|L2,M2> states. Syntax: cg=clebsch_gordan(L,M,L1,M1,L2,M2)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 34-35: Check consistency; implemented by `grumble(L,M,L1,M1,L2,M2)`.
- Lines 37-38: Match the notation to Varshalovich, Section 8.2.1; implemented by `cg=0; c=L; gam=M; a=L1; alp=M1; b=L2; bet=M2`.
- Lines 40-44: Run zero tests (Stage I); implemented by `prefactor_is_nonzero=(a+alp>=0)&&(a-alp>=0)&& (b+bet>=0)&&(b-bet>=0)&& (c+gam>=0)&&(c-gam>=0)&& (gam==alp+bet)`.
- Lines 46-47: Proceed if appropriate; implemented by `if prefactor_is_nonzero`.
- Lines 49-51: Run zero tests (Stage II); implemented by `delta_is_nonzero=(a+b-c>=0)&&(a-b+c>=0)&& (-a+b+c>=0)&&(a+b+c+1>=0)`.
- Lines 53-54: Proceed if appropriate; implemented by `if delta_is_nonzero`.
- Lines 56-57: Run zero tests (Stage III); implemented by `lower_sum_limit=max([alp-a, b+gam-a, 0])`.
- Lines 61-62: Proceed if appropriate; implemented by `if sum_is_nonzero`.
- Lines 64-65: Compile a look-up table of gamma functions; implemented by `gamma(a+b+c+2)=java.math.BigInteger.ONE; gamma(:)=gamma(a+b+c+2)`.
- Lines 70-71: Compute prefactor numerator; implemented by `numer=java.math.BigInteger.valueOf(2*c+1)`.
- Lines 78-79: Compute prefactor denominator; implemented by `denom=gamma(a+b+c+2)`.
- Lines 85-86: Perform floating-point division; implemented by `accur=length(denom.toString)-length(numer.toString)+64`.
- Lines 90-91: Compute the sum; implemented by `z_sum=java.math.BigDecimal.ZERO`.
- Lines 94-95: Compute term numerator; implemented by `numer=java.math.BigInteger.valueOf((-1)^(b+bet+z))`.
- Lines 99-100: Compute term denominator; implemented by `denom=gamma(z+1)`.
- Lines 105-106: Perform floating-point division; implemented by `numer=java.math.BigDecimal(numer); denom=java.math.BigDecimal(denom)`.
- Lines 109-110: Add the term to the total; implemented by `z_sum=z_sum.add(sum_term)`.
- Lines 114-115: Return to double precision; implemented by `cg_sq=z_sum.multiply(z_sum).multiply(prefactor)`.

### Control flow inferred from the code

- Line 47: conditional branch on `prefactor_is_nonzero`.
- Line 54: conditional branch on `delta_is_nonzero`.
- Line 62: conditional branch on `sum_is_nonzero`.
- Line 66: `for` loop over `k=3:(a+b+c+2)`.
- Line 92: `for` loop over `z=lower_sum_limit:upper_sum_limit`.

### Key state/data transformations

- Lines 38: computes `cg` using `cg=0; c=L; gam=M; a=L1; alp=M1; b=L2; bet=M2`.
- Lines 57: computes `lower_sum_limit` using `lower_sum_limit=max([alp-a, b+gam-a, 0])`.
- Lines 58: computes `upper_sum_limit` using `upper_sum_limit=min([c+b+alp, c+b-a, c+gam])`.
- Lines 65: computes `gamma(a+b+c+2)` using `gamma(a+b+c+2)=java.math.BigInteger.ONE; gamma(:)=gamma(a+b+c+2)`.
- Lines 67: computes `gamma(k)` using `gamma(k)=gamma(k-1).multiply(java.math.BigInteger.valueOf(k-1))`.
- Lines 71: computes `numer` using `numer=java.math.BigInteger.valueOf(2*c+1)`.
- Lines 79: computes `denom` using `denom=gamma(a+b+c+2)`.
- Lines 86: computes `accur` using `accur=length(denom.toString)-length(numer.toString)+64`.
- Lines 88: computes `prefactor` using `prefactor=numer.divide(denom,accur,java.math.RoundingMode.HALF_UP)`.
- Lines 91: computes `z_sum` using `z_sum=java.math.BigDecimal.ZERO`.
- Lines 107: computes `sum_term` using `sum_term=numer.divide(denom,accur,java.math.RoundingMode.HALF_UP)`.
- Lines 115: computes `cg_sq` using `cg_sq=z_sum.multiply(z_sum).multiply(prefactor)`.

### Local helper functions

- Line 127: `grumble()` — `function grumble(L,M,L1,M1,L2,M2)`.
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
- a trivial matter for high ranks. This function produces machine
- precision answers up to about L=1e4. A faster implementation for
- low ranks is available in cg_fast.m function.

## Implementation structure

- Clebsch-Gordan coefficient: the coefficient in front of Y(L,M) spheri-
- cal harmonic in the expansion of the product of Y(L1,M1) and Y(L2,M2)
- spherical harmonics. In the more general sense, the coefficient refers
- to the expansion coefficient of |L,M> angular momentum or spin state
- in the product basis of |L1,M1>|L2,M2> states. Syntax:
- cg=clebsch_gordan(L,M,L1,M1,L2,M2)
- L,M,L1,M1,L2,M2 -integer or half-integer indices of
- the angular momentum or spin states
- cg -floating-point (double precision)
- Clebsch-Gordan coefficient
- Note: only some combinations of L,M,L1,M1,L2,M2 are allowed by the pro-
- perties of spherical harmonics and spin states. If inadmissible

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `gamma()`, `double()`.
