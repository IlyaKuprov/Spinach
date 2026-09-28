# examples/relaxation_theory/aniso_diff_test_2.m

- Signature: `aniso_diff_test_2()`

## Purpose

Build and display the Redfield relaxation superoperator for an anisotropically shielded, dipole-coupled `1H`–`13C` pair undergoing anisotropic rotational diffusion. Calculation time: seconds.

## Model and parameters

- The field is set from `2*pi*950.33e6/spin('1H')`; the scalar coupling is `145.0 Hz`.
- The two shielding-tensor principal-value rows are `[10 20 30] ppm` and `[40 50 60] ppm`. Both Euler-angle rows are `[0 pi/4 0]`.
- Rotational-diffusion eigenvalues are `[2.16e8 2.35e8 7.45e8]`, with correlation-time parameter `1./(6*D)`.
- Relaxation is Redfield, with zero equilibrium and `labframe` retention; the basis is `sphten-liouv` with no approximation.

## Calculation

The function creates and bases the Spinach system, evaluates `relaxation(spin_system)`, and prints the full matrix with `disp(full(R))`.
