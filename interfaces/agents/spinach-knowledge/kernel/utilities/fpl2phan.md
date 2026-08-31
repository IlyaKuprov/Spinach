# kernel/utilities/fpl2phan.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/fpl2phan.m`
- Signature: `phan=fpl2phan(rho,coil,dims)`
- Total lines: 59

## Purpose

Returns the image painted within the Fokker-Planck vector by the user-specified spin state. Syntax: phan=fpl2phan(rho,coil,dims)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 25-26: Check consistency; implemented by `grumble(rho,coil,dims)`.
- Lines 28-29: Expose the spin dimension; implemented by `rho=reshape(rho,[numel(coil) prod(dims)])`.
- Lines 31-32: Compute the observable; implemented by `phan=coil'*rho`.
- Lines 34-35: Reshape as needed; implemented by `phan=reshape(phan,dims)`.

### Key state/data transformations

- Lines 29: computes `rho` using `rho=reshape(rho,[numel(coil) prod(dims)])`.
- Lines 32: computes `phan` using `phan=coil'*rho`.

### Local helper functions

- Line 40: `grumble()` — `function grumble(rho,coil,dims)`.
  - Representative operation: `if (~isnumeric(rho))||(size(rho,2)~=1)`.
  - Representative operation: `error('rho must be a column vector.')`.

## Parameters / inputs

- rho -state vector in Fokker-Planck space
- coil -observable state vector in Liouville space
- dims -spatial dimensions of the Fokker-Planck
- problem, a row vector of integers
- Output:
- phan -the image painted by the specified state

## Implementation structure

- Returns the image painted within the Fokker-Planck vector by
- the user-specified spin state. Syntax:
- phan=fpl2phan(rho,coil,dims)
- rho -state vector in Fokker-Planck space
- coil -observable state vector in Liouville space
- dims -spatial dimensions of the Fokker-Planck
- problem, a row vector of integers
- Output:
- phan -the image painted by the specified state
- Check consistency
- Expose the spin dimension
- Compute the observable

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `any()`, `space()`.
