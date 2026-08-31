# kernel/utilities/rspert.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/rspert.m`
- Signature: `[Ep,Vp]=rspert(E0,H1,order)`
- Total lines: 105

## Purpose

Rayleigh-Schrodinger perturbation theory to arbitrary order, Eqs 2.21-2.23 from Stefan Stoll's PhD thesis, with the typo fixed in the numerator of Eq 2.21. Syntax: [Ep,Vp]=rspert(E0,H1,order)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 41-42: Check consistency; implemented by `grumble(E0,H1,order)`.
- Lines 44-45: Auxiliary arrays; implemented by `E=cell(order,1); V=cell(order,1)`.
- Lines 47-48: Reciprocal energy differences; implemented by `Q=1./(E0'-E0); Q(logical(speye(size(Q))))=0`.
- Lines 53-54: First order; implemented by `E{1}=diag(H1); V{1}=Q.*H1`.
- Lines 56-57: Higher orders; implemented by `for k=2:order`.
- Lines 59-60: Expensive; implemented by `R=H1*V{k-1}`.
- Lines 62-63: Eigenvalues; implemented by `E{k}=real(diag(R))`.
- Lines 65-66: Eigenvectors; implemented by `V{k}=R`.
- Lines 74-75: Summation; implemented by `Ep=E0; for n=1:order, Ep=Ep+E{n}; end`.
- Lines 78-79: Normalisation; implemented by `Vp=Vp./sqrt(sum(abs(Vp).^2,1))`.

### Control flow inferred from the code

- Line 49: conditional branch on `any(~isfinite(Q),'all')`.
- Line 57: `for` loop over `k=2:order`.
- Line 67: `for` loop over `m=1:(k-1)`.

### Key state/data transformations

- Lines 45: computes `E` using `E=cell(order,1); V=cell(order,1)`.
- Lines 48: computes `Q` using `Q=1./(E0'-E0); Q(logical(speye(size(Q))))=0`.
- Lines 54: computes `E{1}` using `E{1}=diag(H1); V{1}=Q.*H1`.
- Lines 60: computes `R` using `R=H1*V{k-1}`.
- Lines 63: computes `E{k}` using `E{k}=real(diag(R))`.
- Lines 66: computes `V{k}` using `V{k}=R`.
- Lines 75: computes `Ep` using `Ep=E0; for n=1:order, Ep=Ep+E{n}; end`.
- Lines 76: computes `Vp` using `Vp=eye(size(H1)); for n=1:order, Vp=Vp+V{n}; end`.

### Local helper functions

- Line 84: `grumble()` — `function grumble(E0,H1,order)`.
  - Representative operation: `if (~isnumeric(E0))||(~isreal(E0))||(~iscolumn(E0))`.
  - Representative operation: `error('E0 must be a real column vector.')`.

## Parameters / inputs

- E0 -eigenvalues of H0, a column vector of real
- numbers
- H1 -perturbation, written in the basis that di-
- agonalises H0
- order -order of perturbation theory to be used, 6
- is the sensible maximum

## Outputs

- Ep -eigenvalues of H0+H1 to the specified order,
- a vector of reals, not necessarily sorted in
- the same way as the input
- Vp -normalised eigenvectors of H0+H1 to the spe-
- cified order in perturbation theory, a squa-
- re unitary matrix with eigenvectors in cols
- in the same order as the eigenvalues in Ep
- Notes: there must be no degeneracies in H0; H1 must be Hermitian,
- the theory only converges when norm(H1,2) is much smaller
- than the smallest energy gap in H0; numerical artefacts
- appear beyond sixth order; complexity is linear in the or-
- der and cubic in the matrix dimension.

## Implementation structure

- Rayleigh-Schrodinger perturbation theory to arbitrary order, Eqs
- 2.21-2.23 from Stefan Stoll's PhD thesis, with the typo fixed in
- the numerator of Eq 2.21. Syntax:
- [Ep,Vp]=rspert(E0,H1,order)
- E0 -eigenvalues of H0, a column vector of real
- numbers
- H1 -perturbation, written in the basis that di-
- agonalises H0
- order -order of perturbation theory to be used, 6
- is the sensible maximum
- Ep -eigenvalues of H0+H1 to the specified order,
- a vector of reals, not necessarily sorted in

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `logical()`, `speye()`, `any()`, `iscolumn()`, `ishermitian()`, `isscalar()`.
