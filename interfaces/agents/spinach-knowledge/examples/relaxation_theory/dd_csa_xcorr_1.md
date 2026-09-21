# examples/relaxation_theory/dd_csa_xcorr_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/dd_csa_xcorr_1.m`
- Signature: `dd_csa_xcorr_1()`
- Total lines: 44

## Purpose

Complete Bloch-Redfield-Wangsness relaxation superoperator in a system with two anisotropically shielded nuclei with a dipolar coupling betwe- en them. Spinach relaxation theory module automatically accounts for all cross-correlations (CSA-CSA and DD-CSA cross-correlations are both pre- sent in this case). Dipolar couplings are computed from Cartesian coor- dinates of the two spins. Calculation time: seconds

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Complete Bloch-Redfield-Wangsness relaxation superoperator in a system
- with two anisotropically shielded nuclei with a dipolar coupling betwe-
- en them. Spinach relaxation theory module automatically accounts for all
- cross-correlations (CSA-CSA and DD-CSA cross-correlations are both pre-
- sent in this case). Dipolar couplings are computed from Cartesian coor-
- dinates of the two spins.
- Calculation time: seconds
- Spin system
- Basis set
- Interactions
- Relaxation theory parameters
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `relaxation()`.
