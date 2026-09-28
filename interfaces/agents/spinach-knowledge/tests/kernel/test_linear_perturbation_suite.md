# tests/kernel/test_linear_perturbation_suite.m

- Signature: `result=test_linear_perturbation_suite()`

## Purpose

Tests linear-algebra, angular-momentum, and perturbation utilities. Syntax: result=test_linear_perturbation_suite()

## Physical / mathematical content

- Addition of two spin-`1/2` systems produces one singlet and one triplet. The tests check their multiplicities, projector completeness, and orthonormality.
- A two-level Hermitian system provides reference second-order Rayleigh–Schrödinger and Van Vleck energy shifts. A four-level Hermitian system provides an exact-diagonalisation reference for eighth-order Van Vleck energies.

## Numerical / algorithmic content

- Checks that Rayleigh–Schrödinger eigenvectors are column-normalised and that second- and eighth-order Van Vleck generators are anti-Hermitian.
- Tests analytical Tikhonov inversion with identity data and regularisation matrices: for regularisation parameter `1/2`, the solution is `fit_rhs/(1+reg_param)`. It also checks the reported squared fit error `norm(K*x-y,2)^2` and regularisation signal `norm(D*x,2)^2`.
- Recovers a linear transfer matrix from a full-row-rank set of input/output vector pairs.
- Compares a finite-difference Jacobian with analytical derivatives and checks that its error estimates are finite and non-negative.

## Outputs

- `result` — regression test result with explanatory messages.

## Implementation structure

- Creates a test result for `kernel/linear_perturbation_suite` and records comparisons using `test_close` and `test_true`.
- Exercises `add_spins`, `rspert`, `vvpert`, `tikhoind`, `transfermat`, and `jacobianest` on small reference cases.