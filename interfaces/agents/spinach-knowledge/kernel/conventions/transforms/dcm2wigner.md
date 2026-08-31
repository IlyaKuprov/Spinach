# kernel/conventions/transforms/dcm2wigner.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/dcm2wigner.m`
- Signature: `D=dcm2wigner(dcm)`
- Total lines: 82

## Purpose

Converts a directional cosine matrix into second-rank Wigner function matrix. Syntax: D=dcm2wigner(dcm)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 29-30: Check consistency; implemented by `grumble(dcm)`.
- Lines 32-33: Compute A and B coefficients; implemented by `A=sqrt(0.5*(+dcm(1,1)+1i*dcm(1,2)-1i*dcm(2,1)+dcm(2,2)))`.
- Lines 36-37: Verify amplitudes; implemented by `if abs(A*A'-0.5*(1+dcm(3,3)))+abs(B*B'-0.5*(1-dcm(3,3)))>1e-6`.
- Lines 41-42: Verify phases; implemented by `if abs(A*B+0.5*(dcm(1,3)-1i*dcm(2,3)))+abs(A*B'-0.5*(dcm(3,1)+1i*dcm(3,2)))>1e-6`.
- Lines 49-50: Compute Wigner matrix; implemented by `Z=A*A'-B*B'`.

### Control flow inferred from the code

- Line 37: conditional branch on `abs(A*A'-0.5*(1+dcm(3,3)))+abs(B*B'-0.5*(1-dcm(3,3)))>1e-6`.
- Line 42: conditional branch on `abs(A*B+0.5*(dcm(1,3)-1i*dcm(2,3)))+abs(A*B'-0.5*(dcm(3,1)+1i*dcm(3,2)))>1e-6`.
- Line 45: conditional branch on `abs(A*B+0.5*(dcm(1,3)-1i*dcm(2,3)))+abs(A*B'-0.5*(dcm(3,1)+1i*dcm(3,2)))>1e-6`.

### Key state/data transformations

- Lines 33: computes `A` using `A=sqrt(0.5*(+dcm(1,1)+1i*dcm(1,2)-1i*dcm(2,1)+dcm(2,2)))`.
- Lines 34: computes `B` using `B=sqrt(0.5*(-dcm(1,1)+1i*dcm(1,2)+1i*dcm(2,1)+dcm(2,2)))`.
- Lines 50: computes `Z` using `Z=A*A'-B*B'`.
- Lines 51: computes `D` using `D=[ A^4 2*A^3*B sqrt(6)*A^2*B^2 2*A*B^3 B^4`.

### Local helper functions

- Line 60: `grumble()` — `function grumble(dcm)`.
  - Representative operation: `if (~isnumeric(dcm))||(~isreal(dcm))||(~all(size(dcm)==[3 3]))`.
  - Representative operation: `error('DCM must be a real 3x3 matrix.')`.

## Parameters / inputs

- dcm -directional cosine matrix

## Outputs

- D -matrix of second rank Wigner D functions. Rows
- and columns are sorted by descending ranks:
- [D( 2,2) ... D( 2,-2)
- ... ... ...
- D(-2,2) ... D(-2,-2)]
- Notes: the resulting Wigner matrix is to be used as v=W*v, where v is
- a column vector of irreducible spherical tensor coefficients in
- the following order: T(2,2), T(2,1), T(2,0), T(2,-1), T(2,-2).

## Implementation structure

- Converts a directional cosine matrix into second-rank Wigner function
- matrix. Syntax:
- D=dcm2wigner(dcm)
- dcm -directional cosine matrix
- D -matrix of second rank Wigner D functions. Rows
- and columns are sorted by descending ranks:
- [D( 2,2) ... D( 2,-2)
- ... ... ...
- D(-2,2) ... D(-2,-2)]
- a column vector of irreducible spherical tensor coefficients in
- the following order: T(2,2), T(2,1), T(2,0), T(2,-1), T(2,-2).
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `dcm()`, `all()`.
