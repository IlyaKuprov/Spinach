# tests/kernel/test_dynamic_relaxation_models_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_relaxation_models_suite.m`
- Signature: `result=test_dynamic_relaxation_models_suite()`
- Total lines: 160

## Purpose

Tests dynamic relaxation model helper paths. Syntax: result=test_dynamic_relaxation_models_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Outputs

- result -regression test result with explanatory messages
- The test checks anisotropic and functional T1/T2 rates, damping, Lindblad,
- scalar Redfield, correlation functions, and the serial Redfield include.

## Implementation structure

- Tests dynamic relaxation model helper paths. Syntax:
- result=test_dynamic_relaxation_models_suite()
- result -regression test result with explanatory messages
- The test checks anisotropic and functional T1/T2 rates, damping, Lindblad,
- scalar Redfield, correlation functions, and the serial Redfield include.
- Announce the test target
- State the relaxation-model target of the test
- Check anisotropic tensor rates in the extended T1/T2 model
- Check function-handle rates in the extended T1/T2 model
- Check non-selective damping in Liouville space
- Check one-spin Lindblad relaxation rates
- Check scalar Redfield integral against the closed zero-H0 reference

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_spin_system()`, `rlx_t1_t2()`, `state()`, `test_close()`, `relaxation()`, `unit_state()`, `rlx_scalar()`, `speye()`, `corrfun()`, `nnz()`, `test_true()`.
