# kernel/utilities/wigner.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/wigner.m`
- Signature: `D=wigner(l,alp,bet,gam)`
- Total lines: 107

## Purpose

Wigner D matrices, defined as (Brink & Satchler, Eq 2.13): D=expm(-1i*Lz*alp)*expm(-1i*Ly*bet)*expm(-1i*Lz*gam); where Lx, Ly, Lz are Pauli matrices. Syntax: D=wigner(l,alp,bet,gam)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- l -rank of the Wigner matrix, may be half-integer
- alp -alp Euler angle, radians
- bet -bet Euler angle, radians
- gam -gam Euler angle, radians
- ZYZ convention is used for Euler angles, see Brink and Satchler,
- Figures 1 and 2.

## Outputs

- D -Wigner D matrix with rows and columns sorted
- by descending ranks, for example (l=2):
- [D( 2,2) ... D( 2,-2)
- ... ... ...
- D(-2,2) ... D(-2,-2)]
- The output is to be used as y=D*x, where x is a column vector of
- irreducible spherical tensor coefficients, listed vertically in
- the order: T(2,2), T(2,1), T(2,0), T(2,-1), T(2,-2).

## Implementation structure

- Wigner D matrices, defined as (Brink & Satchler, Eq 2.13):
- D=expm(-1i*Lz*alp)*expm(-1i*Ly*bet)*expm(-1i*Lz*gam);
- where Lx, Ly, Lz are Pauli matrices. Syntax:
- D=wigner(l,alp,bet,gam)
- l -rank of the Wigner matrix, may be half-integer
- alp -alp Euler angle, radians
- bet -bet Euler angle, radians
- gam -gam Euler angle, radians
- ZYZ convention is used for Euler angles, see Brink and Satchler,
- Figures 1 and 2.
- D -Wigner D matrix with rows and columns sorted
- by descending ranks, for example (l=2):

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `pauli()`, `isscalar()`.
