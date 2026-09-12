# kernel/utilities/wigner_fock.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/wigner_fock.m`
- Signature: `W=wigner_fock(rho,alpha)`
- Total lines: 85

## Purpose

Wigner function of a bosonic mode state given as a density matrix in a truncated Fock basis, evaluated at the specified points of the phase space through the displaced parity operator (Royer, Phys. Rev. A 15, 449, 1977): W(alpha)=(2/pi)*trace(rho*D(alpha)*P*D(alpha)') where D(alpha)=expm(alpha*a'-conj(alpha)*a) is the displacement operator and P=expm(1i*pi*a'*a) is the photon number parity ope- rator, both built from

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- Each requested phase-space point costs one `expm` of the displacement generator alpha*a'-conj(alpha)*a in the truncated Fock basis and one trace against the parity operator; there is no averaging, propagation, or quadrature, and the points are evaluated independently in a loop over `numel(alpha)`.
- The only eigenvalue call is the positive-semidefiniteness test in the grumbler; the Hermiticity, unit-trace, and positivity tolerances are sqrt(eps) of the class of `rho`, so single-precision density matrices are accepted.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 46-47: Check consistency; implemented by `grumble(rho,alpha)`.
- Lines 49-50: Annihilation and parity operators in the truncated Fock basis; implemented by `nlevels=size(rho,1); an_op=diag(sqrt(1:(nlevels-1)),1)`.
- Lines 53-54: Displaced parity expectation values at the grid points; implemented by `W=zeros(size(alpha))`.

### Control flow inferred from the code

- Line 55: `for` loop over `n=1:numel(alpha)`.

### Key state/data transformations

- Lines 50: computes `nlevels` using `nlevels=size(rho,1); an_op=diag(sqrt(1:(nlevels-1)),1)`.
- Lines 51: computes `parity` using `parity=diag((-1).^(0:(nlevels-1)))`.
- Lines 54: computes `W` using `W=zeros(size(alpha))`.
- Lines 56: computes `disp_op` using `disp_op=expm(alpha(n)*an_op'-conj(alpha(n))*an_op)`.
- Lines 57: computes `W(n)` using `W(n)=(2/pi)*real(trace(rho*(disp_op*parity*disp_op')))`.

### Local helper functions

- Line 63: `grumble()` — `function grumble(rho,alpha)`.
  - Representative operation: `if (~isnumeric(rho))||(~ismatrix(rho))||(size(rho,1)~=size(rho,2))||(size(rho,1)<2)||(~all(isfinite(rho),'all'))`.
  - Representative operation: `error('rho must be a square matrix of dimension at least 2 with finite elements.')`.

## Parameters / inputs

- rho -density matrix of the mode in the Fock basis
- with the levels in ascending order, [n x n]
- alpha -complex phase space coordinates, an array of
- any size; the real part is the position qua-
- drature and the imaginary part is the momen-
- tum quadrature in the units where a coherent
- state |beta> has W=(2/pi)*exp(-2*|alpha-beta|^2)

## Outputs

- W -Wigner function values, a real array of the
- same size as alpha, normalised to a unit in-
- tegral over the complex plane
- Note: the Fock basis must be large enough for the displaced
- states to fit: with the highest occupied level m of rho,
- n must be well above (sqrt(m)+|alpha|)^2 at every point,
- otherwise the truncated displacement operator distorts
- the values; pad the density matrix with zero rows and
- columns when the state itself lives in a smaller space,
- and check the unit integral on the grid in use.

## Implementation structure

- Wigner function of a bosonic mode state given as a density matrix
- in a truncated Fock basis, evaluated at the specified points of the
- phase space through the displaced parity operator (Royer, Phys.
- Rev. A 15, 449, 1977):
- W(alpha)=(2/pi)*trace(rho*D(alpha)*P*D(alpha)')
- where D(alpha)=expm(alpha*a'-conj(alpha)*a) is the displacement
- operator and P=expm(1i*pi*a'*a) is the photon number parity ope-
- rator, both built from the ladder operators truncated to the di-
- mension of the density matrix. Syntax:
- W=wigner_fock(rho,alpha)
- rho -density matrix of the mode in the Fock basis
- with the levels in ascending order, [n x n]

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `alpha()`, `conj()`, `ismatrix()`, `all()`, `eps()`.
