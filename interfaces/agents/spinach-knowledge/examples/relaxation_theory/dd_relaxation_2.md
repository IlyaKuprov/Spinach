# examples/relaxation_theory/dd_relaxation_2.m

- Signature: `dd_relaxation_2()`

## Purpose

Complete Bloch-Redfield-Wangsness relaxation superoperator in a system with dipolar coupling between spins. The dipolar couplings are computed from Cartesian coordinates of the spins. The result should not depend on the choice of the rotation angles below. Calculation time: seconds

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Complete Bloch-Redfield-Wangsness relaxation superoperator in a system
- with dipolar coupling between spins. The dipolar couplings are computed
- from Cartesian coordinates of the spins. The result should not depend
- on the choice of the rotation angles below.
- Calculation time: seconds
- System specification
- Randomly rotated set of coordinates
- Relaxation theory parameters
- Basis set
- Spinach housekeeping
- Relaxation superoperator
