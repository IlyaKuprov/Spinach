# tests/kernel/test_linear_perturbation_suite.m

## Purpose

Regression test suite for Spinach linear-algebra, angular-momentum, and perturbation-theory utilities. The file checks spin-addition projectors, Rayleigh–Schrödinger and Van Vleck perturbation theory, analytical Tikhonov inversion, transfer matrices, and finite-difference Jacobian estimation against small analytical references.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_linear_perturbation_suite.m>

## Behaviour

The function announces the test target with `fprintf`, initialises a result object via `new_test_result` for the `kernel/linear_perturbation_suite` target, and then runs a sequence of `test_close` and `test_true` checks:

- **Spin addition** — calls `add_spins(1/2,1/2)` and verifies that the returned multiplicities equal `[1 3]` (one singlet and one triplet irrep, tolerance `1e-14`), that the projector completeness relation `projectors{1}*projectors{1}'+projectors{2}*projectors{2}'` equals `eye(4)` (`1e-13`), and that each projector block satisfies `projectors{k}'*projectors{k} = eye(1)` / `eye(3)` orthonormality (`1e-13`).
- **Rayleigh–Schrödinger perturbation theory** — for a two-level system with `base_energy=[0;10]` and an off-diagonal perturbation of strength `0.01`, calls `rspert(base_energy,pert_mat,2)` and checks the second-order energies against the reference `[-pert_strength^2/10; 10+pert_strength^2/10]` (`1e-13`) and that the returned eigenvectors are column-normalised (`1e-14`).
- **Van Vleck perturbation theory** — calls `vvpert(base_energy,pert_mat,2)` on the same two-level system and checks the energies against the same reference (`1e-13`) plus anti-Hermiticity of the generator `vv_gen+vv_gen' = 0` (`1e-14`). A second case uses a 4-level system (`base_energy=[-3;-1;2;5]`) with a complex Hermitian perturbation matrix and order 8; the sorted energies are compared to `sort(real(eig(diag(base_energy)+pert_mat,'vector')))` at tolerance `1e-12`, and the generator anti-Hermiticity is checked at `1e-13`.
- **Analytical Tikhonov inversion** — calls `tikhoind(fit_mat,reg_mat,fit_rhs,reg_param)` with identity data and regularisation matrices, `fit_rhs=[3;6]`, and `reg_param=1/2`; verifies the solution equals `fit_rhs/(1+reg_param)` (`1e-14`), the error signal equals `norm(tikh_ref-fit_rhs,2)^2` (`1e-14`), and the regularisation signal equals `norm(tikh_ref,2)^2` (`1e-14`).
- **Transfer matrix** — builds `amp_inputs=[1 0 1 2;0 1 1 -1]` and `transfer_ref=[2 -1;1/2 3]`, computes `amp_outputs=transfer_ref*amp_inputs`, and checks that `transfermat(amp_inputs,amp_outputs)` recovers `transfer_ref` exactly (`1e-13`).
- **Finite-difference Jacobian** — defines `jac_fun=@(x)[x(1)^2+3*x(2);sin(x(1)*x(2))]` at `jac_point=[2;0.3]`, calls `jacobianest(jac_fun,jac_point)`, and compares the estimated Jacobian to the analytical matrix `[4 3;0.3*cos(0.6) 2*cos(0.6)]` (`1e-6`); additionally asserts via `test_true` that all error estimates are finite and non-negative.

Each check appends an explanatory message to the result object, which is returned to the caller.

## Inputs and outputs

**Syntax**

```matlab
result = test_linear_perturbation_suite()
```

The function takes no inputs.

**Outputs**

- `result` — regression test result object with explanatory messages, accumulated through the `test_close` and `test_true` checks described above.

## References

1. Spinach library test suite, `tests/kernel/test_linear_perturbation_suite.m`, <https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_linear_perturbation_suite.m>
