# kernel/utilities/clebsch_gordan.m

- Signature: `cg=clebsch_gordan(L,M,L1,M1,L2,M2)`

## Purpose

Clebsch-Gordan coefficient: the coefficient in front of Y(L,M) spheri- cal harmonic in the expansion of the product of Y(L1,M1) and Y(L2,M2) spherical harmonics. In the more general sense, the coefficient refers to the expansion coefficient of |L,M> angular momentum or spin state in the product basis of |L1,M1>|L2,M2> states. Syntax: cg=clebsch_gordan(L,M,L1,M1,L2,M2)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

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
