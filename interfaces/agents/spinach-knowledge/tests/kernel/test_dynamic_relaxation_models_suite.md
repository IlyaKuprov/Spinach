# tests/kernel/test_dynamic_relaxation_models_suite.m

- Signature: `result=test_dynamic_relaxation_models_suite()`

## Purpose

Tests dynamic relaxation model helper paths. Syntax: result=test_dynamic_relaxation_models_suite()

## Physical / mathematical content

- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- Builds one-spin Liouville-space examples and checks anisotropic tensor and function-handle T1/T2 rates, damping and Lindblad rates, a zero-H0 scalar Redfield integral, rotational-diffusion correlation-function weights and rates, and the serial Redfield include.
## Outputs

- `result` — regression test result with explanatory messages.
## Implementation structure

- Check anisotropic-tensor and function-handle T1/T2 rates.
- Check unit-state preservation and relaxation rates for damping and one-spin Lindblad models.
- Compare the scalar Redfield integral with its zero-H0 reference; check isotropic, axial, and rhombic correlation-function weights and rates, including species state count.
- Exercise the serial Redfield include and check that its matrix is nonzero and preserves the unit state.