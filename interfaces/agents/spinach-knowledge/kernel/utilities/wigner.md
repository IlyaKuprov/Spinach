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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 41-42: Check consistency; implemented by `grumble(l,alp,bet,gam)`.
- Lines 44-45: Rank 1 hard-coded for speed; implemented by `if l==1`.
- Lines 47-48: Blob components; implemented by `egp1=exp(-1i*1*gam); egm1=exp(-1i*(-1)*gam)`.
- Lines 52-55: This ugly blob runs faster than the generic expm call below; implemented by `D=[eap1*egp1*(1/2+cb/2), -eap1*sb/sqrt(2), eap1*egm1*(1/2-cb/2); egp1*sb/sqrt(2), cb, -egm1*sb/sqrt(2); eam1*egp1*(1/2-cb/2), eam1*sb/sqrt(2), eam1*egm1*(1/2+cb/2)]`.
- Lines 57-58: Rank 2 hard-coded for speed; implemented by `elseif l==2`.
- Lines 60-61: Blob components; implemented by `egp2=exp(-1i*2*gam); egp1=exp(-1i*1*gam)`.
- Lines 67-77: This ugly blob runs faster than the generic expm call below; implemented by `D=[eap2*cbh4*egp2, eap2*(-sb*(1+cb)/2)*egp1, eap2*sqrt(3/8)*sb^2, eap2*(-sb*(1-cb)/2)*egm1, eap2*sbh4*egm2; eap1*(sb*(cb+1)/2)*egp2, eap1*(2*cb-1)*(1+cb)/2*egp1, eap1*(-…`.
- Lines 81-82: Get Pauli matrices; implemented by `L=pauli(2*l+1)`.
- Lines 84-85: Compute Wigner matrix (Brink & Satchler, Eq 2.13); implemented by `D=expm(-1i*L.z*alp)*expm(-1i*L.y*bet)*expm(-1i*L.z*gam)`.

### Control flow inferred from the code

- Line 45: conditional branch on `l==1`.

### Key state/data transformations

- Lines 48: computes `egp1` using `egp1=exp(-1i*1*gam); egm1=exp(-1i*(-1)*gam)`.
- Lines 49: computes `eap1` using `eap1=exp(-1i*1*alp); eam1=exp(-1i*(-1)*alp)`.
- Lines 50: computes `sb` using `sb=sin(bet); cb=cos(bet)`.
- Lines 53-55: computes `D` using `D=[eap1*egp1*(1/2+cb/2), -eap1*sb/sqrt(2), eap1*egm1*(1/2-cb/2); egp1*sb/sqrt(2), cb, -egm1*sb/sqrt(2); eam1*egp1*(1/2-cb/2), eam1*sb/sqrt(2), eam1*egm1*(1/2+cb/2)]`.
- Lines 61: computes `egp2` using `egp2=exp(-1i*2*gam); egp1=exp(-1i*1*gam)`.
- Lines 62: computes `egm1` using `egm1=exp(-1i*(-1)*gam); egm2=exp(-1i*(-2)*gam)`.
- Lines 63: computes `eap2` using `eap2=exp(-1i*2*alp); eap1=exp(-1i*1*alp)`.
- Lines 64: computes `eam1` using `eam1=exp(-1i*(-1)*alp); eam2=exp(-1i*(-2)*alp)`.
- Lines 82: computes `L` using `L=pauli(2*l+1)`.

### Local helper functions

- Line 92: `grumble()` — `function grumble(l,alp,bet,gam)`. Every stink that fights the ventilator thinks it
  - Representative operation: `if (~isnumeric(l))||(~isreal(l))||(~isscalar(l))||(mod(l,1/2)~=0)||(l<0)`.
  - Representative operation: `error('l must be a non-negative real integer or half-integer.')`.

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
