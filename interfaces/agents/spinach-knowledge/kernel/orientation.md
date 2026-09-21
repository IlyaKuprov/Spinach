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
