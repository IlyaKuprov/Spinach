# examples/fundamentals/convention_tests/rotations_2.m

- Signature: `rotations_2()`

## Purpose

Compares two representations of the same rotated spin system: rotate the interaction tensors and coordinates in the input, or keep them fixed and rotate the anisotropic Hamiltonian contribution with Spinach's `orientation` function.

## Method and check

The test generates random chemical-shift and coupling tensors for a 1H–15N pair at 14.1 T, using the `sphten-liouv` basis with no approximation. In the kernel-level construction, it evaluates the Hamiltonian at Euler angles [1,2,3]. In the input-level construction, it rotates both shift tensors, the coupling tensor, and the coordinates by the corresponding DCM, then evaluates the Hamiltonian at zero orientation. The 1-norm of the Hamiltonian difference must be no greater than 10⁻³.
