# examples/relaxation_theory/aniso_diff_test_1.m

- Signature: `aniso_diff_test_1()`

## Purpose

Build and display a Redfield relaxation superoperator for a dipole-coupled `1H`–`13C` pair with anisotropic rotational diffusion and zero chemical-shift anisotropy. Calculation time: seconds.

## Model and parameters

- The field is set from `2*pi*950.33e6/spin('1H')`.
- Both shielding tensors have zero principal values and zero Euler angles.
- The coordinates (Å) are `[1.13 0 0]*R` and the origin, where `R=euler2dcm(1,2,3)`; the scalar coupling is `145.0 Hz`.
- Rotational-diffusion eigenvalues are `[2.16e8 2.35e8 7.45e8]`; the correlation-time parameter is `1./(6*D)`.
- Relaxation is Redfield, with zero equilibrium and `labframe` retention. The basis is `sphten-liouv` with no approximation.

## Calculation

The function creates and bases the Spinach system, evaluates `relaxation(spin_system)`, and prints the full matrix with `disp(full(R))`.
