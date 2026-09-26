# examples/relaxation_theory/from_md/ngce_test.m

- Signature: `ngce_test()`

## Purpose

Test of the numerical integral route to the Redfield relaxation superoperator against the analytical results for isotropic rota- tional diffusion. Calculation time: minutes.

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- Test of the numerical integral route to the Redfield relaxation
- superoperator against the analytical results for isotropic rota-
- tional diffusion.
- Calculation time: minutes.
- System specification
- Relaxation theory parameters
- Basis set
- Spinach housekeeping
- Get Redfield relaxation matrix
- Get lab frame Hamiltonian components
- Get a random walk on a sphere
- Get Hamiltonian trajectory
