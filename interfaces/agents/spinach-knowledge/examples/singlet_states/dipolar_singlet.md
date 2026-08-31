# examples/singlet_states/dipolar_singlet.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/singlet_states/dipolar_singlet.m`
- Signature: `dipolar_singlet()`
- Total lines: 48

## Purpose

A demonstration that the two-spin singet state is immune to dipolar relaxation. Full Redfield superoperator for dipolar relaxation in liquid state is computed and the norm of its action on a singlet state is printed to the console. Calculation time: seconds

## Physical / mathematical content

- Long-lived singlet-state examples. The central concept is symmetry-protected or nearly symmetry-protected two-spin order that relaxes much more slowly than ordinary Zeeman magnetisation. Files here often analyse singlet-triplet subspaces, state conversion sequences, and relaxation leakage channels.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: System specification; implemented by `sys.magnet=14.1`.
- Lines 19-20: Relaxation theory parameters; implemented by `inter.relaxation={'redfield'}`.
- Lines 25-26: Relaxation superoperator accuracy; implemented by `sys.tols.rlx_integration=1e-5`.
- Lines 29-30: Proximity cut-off; implemented by `sys.tols.prox_cutoff=4.0`.
- Lines 32-33: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 36-37: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 40-41: Relaxation superoperator; implemented by `R=relaxation(spin_system)`.
- Lines 43-44: Action on a singlet state; implemented by `S=singlet(spin_system,1,2); S=S/norm(S)`.

### Key state/data transformations

- Lines 14: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 15: computes `sys.isotopes` using `sys.isotopes={'1H','1H'}`.
- Lines 16: computes `inter.coordinates` using `inter.coordinates={[0.0 0.0 0.0]`.
- Lines 20: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 21: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 22: computes `inter.rlx_keep` using `inter.rlx_keep='labframe'`.
- Lines 23: computes `inter.tau_c` using `inter.tau_c={5e-9}`.
- Lines 26: computes `sys.tols.rlx_integration` using `sys.tols.rlx_integration=1e-5`.
- Lines 27: computes `sys.tols.rlx_zero` using `sys.tols.rlx_zero=1e-5`.
- Lines 30: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 33: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 34: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 37: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 41: computes `R` using `R=relaxation(spin_system)`.
- Lines 44: computes `S` using `S=singlet(spin_system,1,2); S=S/norm(S)`.

## Implementation structure

- A demonstration that the two-spin singet state is immune
- to dipolar relaxation. Full Redfield superoperator for
- dipolar relaxation in liquid state is computed and the
- norm of its action on a singlet state is printed to the
- console.
- Calculation time: seconds
- System specification
- Relaxation theory parameters
- Relaxation superoperator accuracy
- Proximity cut-off
- Basis set
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `relaxation()`, `singlet()`, `num2str()`.
