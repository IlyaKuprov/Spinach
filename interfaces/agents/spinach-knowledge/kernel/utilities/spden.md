# kernel/utilities/spden.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/spden.m`
- Signature: `J=spden(L,D,omega)`
- Total lines: 62

## Purpose

Lorentzian spectral density function for rotational diffusion at the user-specified frequency. Syntax: J=spden(L,D,omega)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 25-26: Check consistency; implemented by `grumble(L,D,omega)`.
- Lines 28-29: Get the correlation time; implemented by `tau_c=1/(L*(L+1)*D)`.
- Lines 31-32: Get the spectral density; implemented by `J=(tau_c/(2*L+1))/(1+(tau_c*omega)^2)`.

### Key state/data transformations

- Lines 29: computes `tau_c` using `tau_c=1/(L*(L+1)*D)`.
- Lines 32: computes `J` using `J=(tau_c/(2*L+1))/(1+(tau_c*omega)^2)`.

### Local helper functions

- Line 37: `grumble()` — `function grumble(L,D,omega)`.
  - Representative operation: `if (~isnumeric(L))||(~isscalar(L))|| (~isreal(L))||(L<1)||(mod(L,1)~=0)`.
  - Representative operation: `(~isreal(L))||(L<1)||(mod(L,1)~=0)`.

## Parameters / inputs

- L -spherical rank, use 2 for common NMR
- mechanisms such as dipolar relaxation
- D -rotational diffusion coefficient, s^{-1}
- omega -frequency, rad/s

## Outputs

- J -spectral density function value

## Implementation structure

- Lorentzian spectral density function for rotational
- diffusion at the user-specified frequency. Syntax:
- J=spden(L,D,omega)
- L -spherical rank, use 2 for common NMR
- mechanisms such as dipolar relaxation
- D -rotational diffusion coefficient, s^{-1}
- omega -frequency, rad/s
- J -spectral density function value
- Check consistency
- Get the correlation time
- Get the spectral density
- Consistency enforcement

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isscalar()`.
