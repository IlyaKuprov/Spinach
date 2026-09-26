# kernel/orientation.m

- Signature: `H=orientation(Q,euler_angles)`

## Purpose

Anisotropic part of the Hamiltonian for a specific spin system orientation. Syntax: H=orientation(Q,euler_angles)

## Physical / mathematical content

## Numerical / algorithmic content

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
