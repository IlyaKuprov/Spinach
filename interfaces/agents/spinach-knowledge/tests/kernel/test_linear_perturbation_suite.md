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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Linear algebra and perturbation utilities\n')`.
- Lines 20-23: State the utility target of the test; implemented by `result=new_test_result('kernel/linear_perturbation_suite', 'Linear algebra and perturbation utilities', 'Small analytical linear-algebra cases must match exact matrix an…`.
- Lines 25-26: Check spin-half addition into singlet and triplet irreducible blocks; implemented by `[multiplicities,projectors]=add_spins(1/2,1/2)`.
- Lines 41-42: Check second-order perturbation energy shifts for a two-level system; implemented by `base_energy=[0;10]`.
- Lines 54-55: Check Van Vleck perturbation theory on the same two-level system; implemented by `[vv_energy,vv_gen]=vvpert(base_energy,pert_mat,2)`.
- Lines 63-64: Check higher-order Van Vleck perturbation theory against diagonalisation; implemented by `base_energy=[-3;-1;2;5]`.
- Lines 78-79: Check the analytical indefinite Tikhonov solution for identity operators; implemented by `fit_mat=eye(2)`.
- Lines 95-96: Check transfer matrix recovery from more vector pairs than dimensions; implemented by `amp_inputs=[1 0 1 2;0 1 1 -1]`.
- Lines 104-105: Check finite-difference Jacobian estimation against analytical derivatives; implemented by `jac_fun=@(x)[x(1)^2+3*x(2);sin(x(1)*x(2))]`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/linear_perturbation_suite', 'Linear algebra and perturbation utilities', 'Small analytical linear-algebra cases must match exact matrix an…`.
- Lines 26: computes `[multiplicities,projectors]` using `[multiplicities,projectors]=add_spins(1/2,1/2)`.
- Lines 30: computes `projector_sum` using `projector_sum=projectors{1}*projectors{1}'+projectors{2}*projectors{2}'`.
- Lines 42: computes `base_energy` using `base_energy=[0;10]`.
- Lines 43: computes `pert_strength` using `pert_strength=0.01`.
- Lines 44: computes `pert_mat` using `pert_mat=[0 pert_strength;pert_strength 0]`.
- Lines 45: computes `pert_ref` using `pert_ref=[-pert_strength^2/10;10+pert_strength^2/10]`.
- Lines 46: computes `[rs_energy,rs_vectors]` using `[rs_energy,rs_vectors]=rspert(base_energy,pert_mat,2)`.
- Lines 55: computes `[vv_energy,vv_gen]` using `[vv_energy,vv_gen]=vvpert(base_energy,pert_mat,2)`.
- Lines 70: computes `exact_energy` using `exact_energy=sort(real(eig(diag(base_energy)+pert_mat,'vector')))`.
- Lines 79: computes `fit_mat` using `fit_mat=eye(2)`.
- Lines 80: computes `reg_mat` using `reg_mat=eye(2)`.
- Lines 81: computes `fit_rhs` using `fit_rhs=[3;6]`.
- Lines 82: computes `reg_param` using `reg_param=1/2`.
- Lines 83: computes `[tikh_x,tikh_err,tikh_reg]` using `[tikh_x,tikh_err,tikh_reg]=tikhoind(fit_mat,reg_mat,fit_rhs,reg_param)`.
- Lines 84: computes `tikh_ref` using `tikh_ref=fit_rhs/(1+reg_param)`.
- Lines 96: computes `amp_inputs` using `amp_inputs=[1 0 1 2;0 1 1 -1]`.
- Lines 97: computes `transfer_ref` using `transfer_ref=[2 -1;1/2 3]`.

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
