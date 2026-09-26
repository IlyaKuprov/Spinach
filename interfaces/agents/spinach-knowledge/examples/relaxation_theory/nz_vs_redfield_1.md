# examples/relaxation_theory/nz_vs_redfield_1.m

- Signature: `nz_vs_redfield_1()`

## Purpose

Nakajima-Zwanzig relaxation theory against Redfield theory for a two-spin system with dipolar and CSA cross-correlations. The three superoperators compared are the off-shell NZ kernel (resolvent form), the on-shell NZ kernel (back-rotated form), and Redfield theory. The on-shell kernel at zero shift reproduces Redfield theory exactly; the off-shell kernel agrees with Redfield theory on the zero-frequency subspace of 

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Nakajima-Zwanzig relaxation theory against Redfield theory for a
- two-spin system with dipolar and CSA cross-correlations. The three
- superoperators compared are the off-shell NZ kernel (resolvent form),
- the on-shell NZ kernel (back-rotated form), and Redfield theory. The
- on-shell kernel at zero shift reproduces Redfield theory exactly; the
- off-shell kernel agrees with Redfield theory on the zero-frequency
- subspace of the coherent Liouvillian and differs in first order in
- omega*tau_c on coherences; a lifetime shift suppresses all rates by
- pushing the kernel off the real axis.
- Calculation time: minutes
- Magnet and isotopes
- Chemical shielding tensors
