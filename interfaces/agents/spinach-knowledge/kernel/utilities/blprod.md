# kernel/utilities/blprod.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/blprod.m`
- Signature: `[X1_AB,X2_AB]=blprod(A,B)`
- Total lines: 66

## Purpose

Extension of Blicharski's tensor invariants into scalar products of different spin interaction tensors using polarisation identi- ties. Syntax: [X1_AB,X2_AB]=blprod(A,B)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 28-29: Check consistency; implemented by `grumble(A,B)`.
- Lines 31-32: Components; implemented by `[LsqAmB,DsqAmB]=blinv(A-B)`.
- Lines 35-36: Polarisation identity, first rank; implemented by `X1_AB=(LsqApB-LsqAmB)/4`.
- Lines 38-39: Polarisation identity, second rank; implemented by `X2_AB=(DsqApB-DsqAmB)/4`.

### Key state/data transformations

- Lines 32: computes `[LsqAmB,DsqAmB]` using `[LsqAmB,DsqAmB]=blinv(A-B)`.
- Lines 33: computes `[LsqApB,DsqApB]` using `[LsqApB,DsqApB]=blinv(A+B)`.
- Lines 36: computes `X1_AB` using `X1_AB=(LsqApB-LsqAmB)/4`.
- Lines 39: computes `X2_AB` using `X2_AB=(DsqApB-DsqAmB)/4`.

### Local helper functions

- Line 44: `grumble()` — `function grumble(A,B)`. This is probably the greatest irony ever: we fought a
  - Representative operation: `if (~isnumeric(A))||(~isreal(A))|| (~ismatrix(A))||any(size(A)~=[3 3])`.
  - Representative operation: `(~ismatrix(A))||any(size(A)~=[3 3])`.

## Parameters / inputs

- A -a real 3x3 matrix
- B -a real 3x3 matrix

## Outputs

- X1_AB -cross-correlation amplitude, first rank
- X2_AB -cross-correlation amplitude, second rank
- Note: this function is not sensitive to the isotropic components
- of A and B tensors.

## Implementation structure

- Extension of Blicharski's tensor invariants into scalar products
- of different spin interaction tensors using polarisation identi-
- ties. Syntax:
- [X1_AB,X2_AB]=blprod(A,B)
- A - a real 3x3 matrix
- B - a real 3x3 matrix
- X1_AB - cross-correlation amplitude, first rank
- X2_AB - cross-correlation amplitude, second rank
- Note: this function is not sensitive to the isotropic components
- of A and B tensors.
- Check consistency
- Components

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `blinv()`, `ismatrix()`, `any()`.
