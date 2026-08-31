# kernel/utilities/xyz2hfc.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/xyz2hfc.m`
- Signature: `A=xyz2hfc(exyz,nxyz,isotope)`
- Total lines: 80

## Purpose

Converts point electron and nuclear coordinates into a hyper- fine interaction tensor. Syntax: A=xyz2hfc(exyz,nxyz,isotope)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 35-36: Check consistency; implemented by `grumble(exyz,nxyz,isotope)`.
- Lines 38-39: Fundamental constants; implemented by `hbar=1.054571730e-34`.
- Lines 42-43: Get magnetogyric ratios; implemented by `gamma_n=spin(isotope)`.
- Lines 45-46: Set the origin; implemented by `nxyz=nxyz-exyz`.
- Lines 48-49: Collect fundamental constants; implemented by `C=10^4*gamma_n*hbar*mu0/(4*pi*(1e-10)^3)`.
- Lines 51-52: Compute the dipolar matrix; implemented by `D=3*(nxyz'*nxyz)/norm(nxyz,2)^5-eye(3)/norm(nxyz,2)^3`.
- Lines 54-55: Compute the dipolar coupling matrix; implemented by `A=C*D`.

### Key state/data transformations

- Lines 39: computes `hbar` using `hbar=1.054571730e-34`.
- Lines 40: computes `mu0` using `mu0=4*pi*1e-7`.
- Lines 43: computes `gamma_n` using `gamma_n=spin(isotope)`.
- Lines 46: computes `nxyz` using `nxyz=nxyz-exyz`.
- Lines 49: computes `C` using `C=10^4*gamma_n*hbar*mu0/(4*pi*(1e-10)^3)`.
- Lines 52: computes `D` using `D=3*(nxyz'*nxyz)/norm(nxyz,2)^5-eye(3)/norm(nxyz,2)^3`.
- Lines 55: computes `A` using `A=C*D`.

### Local helper functions

- Line 60: `grumble()` — `function grumble(exyz,nxyz,isotope)`.
  - Representative operation: `if (~isnumeric(exyz))||(~isreal(exyz))||(~isequal(size(exyz),[1 3]))`.
  - Representative operation: `error('exyz must be a 1x3 real row vector.')`.

## Parameters / inputs

- exyz -Cartesian coordinates of the electron,
- a 1x3 row vector in Angstrom
- nxyz -Cartesian coordinates of the nucleus,
- a 1x3 row vector in Angstrom
- isotope -isitope specification, e.g. '13C'

## Outputs

- A -hyperfine coupling tensor, Gauss
- Note: Gauss units are used for hyperfine couplings because
- they do not depend on the electron g-tensor.
- Note: the tensor returned is the one that enters the spin
- Hamiltonian as S*A*I; it does not scale with the num-
- ber of unpaired electrons because the electron spin
- operator already carries that magnitude.

## Implementation structure

- Converts point electron and nuclear coordinates into a hyper-
- fine interaction tensor. Syntax:
- A=xyz2hfc(exyz,nxyz,isotope)
- exyz -Cartesian coordinates of the electron,
- a 1x3 row vector in Angstrom
- nxyz -Cartesian coordinates of the nucleus,
- isotope -isitope specification, e.g. '13C'
- A -hyperfine coupling tensor, Gauss
- Note: Gauss units are used for hyperfine couplings because
- they do not depend on the electron g-tensor.
- Note: the tensor returned is the one that enters the spin
- Hamiltonian as S*A*I; it does not scale with the num-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spin()`, `isequal()`, `ischar()`.
