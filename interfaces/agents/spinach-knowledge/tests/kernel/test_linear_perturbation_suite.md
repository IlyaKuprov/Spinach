# tests/kernel/test_linear_perturbation_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_linear_perturbation_suite.m`
- Signature: `result=test_linear_perturbation_suite()`
- Total lines: 117

## Purpose

Tests linear-algebra, angular-momentum, and perturbation utilities. Syntax: result=test_linear_perturbation_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.
- The code contains an inverse-problem or ill-conditioning aspect and therefore introduces explicit regularisation, model selection, or stabilisation logic.

## Outputs

- result -regression test result with explanatory messages
- The test checks spin-addition projectors, Rayleigh-Schrödinger
- and Van Vleck perturbation theory, analytical Tikhonov inversion,
- transfer matrices, and finite-difference Jacobian estimation.

## Implementation structure

- Tests linear-algebra, angular-momentum, and perturbation utilities. Syntax:
- result=test_linear_perturbation_suite()
- result -regression test result with explanatory messages
- The test checks spin-addition projectors, Rayleigh-Schrödinger
- and Van Vleck perturbation theory, analytical Tikhonov inversion,
- transfer matrices, and finite-difference Jacobian estimation.
- Announce the test target
- State the utility target of the test
- Check spin-half addition into singlet and triplet irreducible blocks
- Check second-order perturbation energy shifts for a two-level system
- Check Van Vleck perturbation theory on the same two-level system
- Check higher-order Van Vleck perturbation theory against diagonalisation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `add_spins()`, `test_close()`, `rspert()`, `vvpert()`, `tikhoind()`, `transfermat()`, `jacobianest()`, `test_true()`, `all()`, `jac_err()`.
