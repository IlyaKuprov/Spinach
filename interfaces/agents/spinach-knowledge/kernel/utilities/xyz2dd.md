# kernel/utilities/xyz2dd.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/xyz2dd.m`
- Signature: `[d,alp,bet,gam,M]=xyz2dd(r1,r2,isotope1,isotope2)`
- Total lines: 97

## Purpose

Converts coordinate specification of the dipolar interaction into the dipolar interaction constant, three Euler angles, and the dipolar interaction matrix. Syntax: [d,alp,bet,gam,m]=xyz2dd(r1,r2,isotope1,isotope2)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 39-40: Check consistency; implemented by `grumble(r1,r2,isotope1,isotope2)`.
- Lines 42-43: Fundamental constants; implemented by `hbar=1.054571628e-34; mu0=4*pi*1e-7`.
- Lines 45-46: Get the distance; implemented by `distance=norm(r2-r1,2)`.
- Lines 48-49: Get the ort; implemented by `ort=(r2-r1)/distance`.
- Lines 51-52: Get the dipolar interaction constant; implemented by `d=spin(isotope1)*spin(isotope2)*hbar*mu0/(4*pi*(distance*1e-10)^3)`.
- Lines 54-55: Get the Euler angles; implemented by `[alp,bet,~]=cart2sph(ort(1),ort(2),ort(3)); bet=pi/2-bet; gam=0`.
- Lines 57-58: Compute the matrix if needed; implemented by `if nargout>4`.
- Lines 60-61: Get the dipolar coupling matrix; implemented by `M=d*[1-3*ort(1)*ort(1) -3*ort(1)*ort(2) -3*ort(1)*ort(3)`.
- Lines 65-66: Clean up rounding errors; implemented by `M=(M+M')/2; M=M-eye(3)*trace(M)/3`.

### Control flow inferred from the code

- Line 58: conditional branch on `nargout>4`.

### Key state/data transformations

- Lines 43: computes `hbar` using `hbar=1.054571628e-34; mu0=4*pi*1e-7`.
- Lines 46: computes `distance` using `distance=norm(r2-r1,2)`.
- Lines 49: computes `ort` using `ort=(r2-r1)/distance`.
- Lines 52: computes `d` using `d=spin(isotope1)*spin(isotope2)*hbar*mu0/(4*pi*(distance*1e-10)^3)`.
- Lines 55: computes `[alp,bet,~]` using `[alp,bet,~]=cart2sph(ort(1),ort(2),ort(3)); bet=pi/2-bet; gam=0`.
- Lines 61: computes `M` using `M=d*[1-3*ort(1)*ort(1) -3*ort(1)*ort(2) -3*ort(1)*ort(3)`.

### Local helper functions

- Line 73: `grumble()` — `function grumble(r1,r2,isotope1,isotope2)`.
  - Representative operation: `if (~isnumeric(r1))||(~isreal(r1))||(numel(r1)~=3)`.
  - Representative operation: `error('r1 must be a three-element real vector.')`.

## Parameters / inputs

- r1,2 -3-element vectors of spin coordinates
- in Angstroms
- isotope1,2 -isotope specification strings, e.g.
- '13C'.

## Outputs

- d -dipolar coupling constant, rad/s
- alp -alpha Euler angle, radians
- bet -beta Euler angle, radians
- gam -gamma Euler angle, radians
- M -dipolar interaction tensor, rad/s
- N.B. Euler angles are not uniquely defined for the orientati-
- on of axial interactions (gamma angle can be anything).
- N.B. free-particle magnetogyric ratios are used, use xyz2hfc.m
- if you have electrons in the system.

## Implementation structure

- Converts coordinate specification of the dipolar interaction
- into the dipolar interaction constant, three Euler angles,
- and the dipolar interaction matrix. Syntax:
- [d,alp,bet,gam,m]=xyz2dd(r1,r2,isotope1,isotope2)
- r1,2 -3-element vectors of spin coordinates
- in Angstroms
- isotope1,2 -isotope specification strings, e.g.
- '13C'.
- d -dipolar coupling constant, rad/s
- alp -alpha Euler angle, radians
- bet -beta Euler angle, radians
- gam -gamma Euler angle, radians

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spin()`, `cart2sph()`, `ort()`, `all()`, `ischar()`.
