# kernel/utilities/blinv.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/blinv.m`
- Signature: `[Lsq,Dsq]=blinv(A)`
- Total lines: 59

## Purpose

Blicharski's relaxation theory invariants, as given by Equations 20-21 in http://doi.org/10.1515/zna-1972-1012, with an error and a typo corrected in Equation 21. Syntax: [Lsq,Dsq]=blinv(A) where A is the interaction matrix. This function is not sensitive to the trace of the matrix. Parameters: A -a real 3x3 matrix

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 24-25: Check consistency; implemented by `grumble(A)`.
- Lines 27-30: First rank invariant; implemented by `Lsq=(A(1,2)-A(2,1))^2+ (A(1,3)-A(3,1))^2+ (A(2,3)-A(3,2))^2`.
- Lines 32-41: Second rank invariant; implemented by `Dsq=A(1,1)^2+ A(2,2)^2+ A(3,3)^2- A(1,1)*A(2,2)- A(1,1)*A(3,3)- A(2,2)*A(3,3)+ (3/4)*((A(1,2)+A(2,1))^2+ (A(1,3)+A(3,1))^2+ (A(2,3)+A(3,2))^2)`.

### Key state/data transformations

- Lines 28-30: computes `Lsq` using `Lsq=(A(1,2)-A(2,1))^2+ (A(1,3)-A(3,1))^2+ (A(2,3)-A(3,2))^2`.
- Lines 33-41: computes `Dsq` using `Dsq=A(1,1)^2+ A(2,2)^2+ A(3,3)^2- A(1,1)*A(2,2)- A(1,1)*A(3,3)- A(2,2)*A(3,3)+ (3/4)*((A(1,2)+A(2,1))^2+ (A(1,3)+A(3,1))^2+ (A(2,3)+A(3,2))^2)`.

### Local helper functions

- Line 46: `grumble()` — `function grumble(A)`. His fight came to an abrupt end with the realisation that truth is ir- relevant, and all that matters in society is belief. Those that see the
  - Representative operation: `if (~isnumeric(A))||(~isreal(A))|| (~ismatrix(A))||any(size(A)~=[3 3])`.
  - Representative operation: `(~ismatrix(A))||any(size(A)~=[3 3])`.

## Outputs

- Lsq -first rank invariant
- Dsq -second rank invariant

## Implementation structure

- Blicharski's relaxation theory invariants, as given by Equations
- 20-21 in http://doi.org/10.1515/zna-1972-1012, with an error and
- a typo corrected in Equation 21. Syntax:
- [Lsq,Dsq]=blinv(A)
- where A is the interaction matrix. This function is not sensitive
- to the trace of the matrix. Parameters:
- A - a real 3x3 matrix
- Lsq - first rank invariant
- Dsq - second rank invariant
- Check consistency
- First rank invariant
- Second rank invariant

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ismatrix()`, `any()`.
