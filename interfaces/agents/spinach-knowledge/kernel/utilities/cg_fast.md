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
