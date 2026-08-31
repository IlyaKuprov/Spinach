# kernel/orientation.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/orientation.m`
- Signature: `H=orientation(Q,euler_angles)`
- Total lines: 80

## Purpose

Anisotropic part of the Hamiltonian for a specific spin system orientation. Syntax: H=orientation(Q,euler_angles)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 31-32: Check consistency; implemented by `grumble(Q,euler_angles)`.
- Lines 34-35: Get started; implemented by `H=sparse(0)`.
- Lines 37-38: Loop over the spherical ranks; implemented by `for r=1:numel(Q)`.
- Lines 40-43: Compute the Wigner matrix; implemented by `D=wigner(r,euler_angles(1), euler_angles(2), euler_angles(3))`.
- Lines 45-46: Update the Hamiltonian; implemented by `for k=1:(2*r+1)`.
- Lines 56-57: Clean up the Hamiltonian; implemented by `H=(H+H')/2`.

### Control flow inferred from the code

- Line 38: `for` loop over `r=1:numel(Q)`.
- Line 46: `for` loop over `k=1:(2*r+1)`.
- Line 47: `for` loop over `m=1:(2*r+1)`.
- Line 48: conditional branch on `nnz(Q{r}{k,m})>0`.

### Key state/data transformations

- Lines 35: computes `H` using `H=sparse(0)`.
- Lines 41-43: computes `D` using `D=wigner(r,euler_angles(1), euler_angles(2), euler_angles(3))`.

### Local helper functions

- Line 62: `grumble()` — `function grumble(Q,euler_angles)`. "Can you explain how Spinach calculates protein structures?"
  - Representative operation: `if ~iscell(Q)`.
  - Representative operation: `error('Q parameter must be a cell array.')`.

## Parameters / inputs

- Q -rotational basis as returned by
- hamiltonian.m function
- euler_angles -a 1x3 vector specifying Euler
- angles (radians) relative to the
- input orientation
- Output:
- H -anisotropic part of the Hamiltonian
- for the specified Euler angles
- Note: this function may be used in both Hilbert and Liouville
- space because the H -> [H, ] adjoint map is linear.
- TODO: efficient sparse summation.

## Implementation structure

- Anisotropic part of the Hamiltonian for a specific spin system
- orientation. Syntax:
- H=orientation(Q,euler_angles)
- Q - rotational basis as returned by
- hamiltonian.m function
- euler_angles - a 1x3 vector specifying Euler
- angles (radians) relative to the
- input orientation
- Output:
- H - anisotropic part of the Hamiltonian
- for the specified Euler angles
- Note: this function may be used in both Hilbert and Liouville

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `wigner()`, `euler_angles()`, `nnz()`, `iscell()`.
