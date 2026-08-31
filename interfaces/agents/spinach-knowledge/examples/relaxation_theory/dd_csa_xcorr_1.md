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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: Spin system; implemented by `sys.isotopes={'1H','13C'}`.
- Lines 17-18: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 21-22: Interactions; implemented by `sys.magnet=14.1`.
- Lines 30-31: Relaxation theory parameters; implemented by `inter.relaxation={'redfield'}`.
- Lines 36-37: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 40-41: Relaxation superoperator; implemented by `disp(full(relaxation(spin_system)))`.

### Key state/data transformations

- Lines 15: computes `sys.isotopes` using `sys.isotopes={'1H','13C'}`.
- Lines 18: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 19: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 22: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 23: computes `inter.zeeman.eigs` using `inter.zeeman.eigs={[7 15 -22]`.
- Lines 25: computes `inter.zeeman.euler` using `inter.zeeman.euler={[pi/3 pi/4 pi/5]`.
- Lines 27: computes `inter.coordinates` using `inter.coordinates={[0.00 0.00 0.00]`.
- Lines 31: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 32: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 33: computes `inter.rlx_keep` using `inter.rlx_keep='labframe'`.
- Lines 34: computes `inter.tau_c` using `inter.tau_c={1e-9}`.
- Lines 37: computes `spin_system` using `spin_system=create(sys,inter)`.

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
