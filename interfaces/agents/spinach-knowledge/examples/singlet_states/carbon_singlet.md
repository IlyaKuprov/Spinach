# examples/singlet_states/carbon_singlet.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/singlet_states/carbon_singlet.m`
- Signature: `carbon_singlet()`
- Total lines: 53

## Purpose

Singlet relaxation rate for the two triple bond carbons in cis-dimethylbut-2-ynedioate. Magnetic parameters com- puted with DFT. Calculation time: seconds

## Physical / mathematical content

- Long-lived singlet-state examples. The central concept is symmetry-protected or nearly symmetry-protected two-spin order that relaxes much more slowly than ordinary Zeeman magnetisation. Files here often analyse singlet-triplet subspaces, state conversion sequences, and relaxation leakage channels.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: System specification; implemented by `sys.magnet=14.1`.
- Lines 23-24: Relaxation theory parameters; implemented by `inter.relaxation={'redfield'}`.
- Lines 29-30: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 33-34: Relaxation superoperator accuracy; implemented by `sys.tols.rlx_integration=1e-5`.
- Lines 37-38: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 41-42: Relaxation superoperator; implemented by `R=relaxation(spin_system)`.
- Lines 44-45: Action on longitudinal magnetization; implemented by `Sz=state(spin_system,{'Lz'},{1}); Sz=Sz/norm(Sz)`.
- Lines 48-49: Action on a singlet state; implemented by `S=singlet(spin_system,1,2); S=S/norm(S)`.

### Key state/data transformations

- Lines 12: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 13: computes `sys.isotopes` using `sys.isotopes={'13C','13C'}`.
- Lines 14: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=[29.13 0.00 0.00`.
- Lines 17: computes `inter.zeeman.matrix{2}` using `inter.zeeman.matrix{2}=[29.13 0.00 0.00`.
- Lines 20: computes `inter.coordinates` using `inter.coordinates={[0.000 0.609 0.298]`.
- Lines 24: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 25: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 26: computes `inter.rlx_keep` using `inter.rlx_keep='labframe'`.
- Lines 27: computes `inter.tau_c` using `inter.tau_c={100e-12}`.
- Lines 30: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 31: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 34: computes `sys.tols.rlx_integration` using `sys.tols.rlx_integration=1e-5`.
- Lines 35: computes `sys.tols.rlx_zero` using `sys.tols.rlx_zero=1e-5`.
- Lines 38: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 42: computes `R` using `R=relaxation(spin_system)`.
- Lines 45: computes `Sz` using `Sz=state(spin_system,{'Lz'},{1}); Sz=Sz/norm(Sz)`.
- Lines 49: computes `S` using `S=singlet(spin_system,1,2); S=S/norm(S)`.

## Implementation structure

- Singlet relaxation rate for the two triple bond carbons
- in cis-dimethylbut-2-ynedioate. Magnetic parameters com-
- puted with DFT.
- Calculation time: seconds
- System specification
- Relaxation theory parameters
- Basis set
- Relaxation superoperator accuracy
- Spinach housekeeping
- Relaxation superoperator
- Action on longitudinal magnetization
- Action on a singlet state

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `relaxation()`, `state()`, `report()`, `num2str()`, `singlet()`.
