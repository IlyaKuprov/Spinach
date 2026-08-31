# examples/singlet_states/decoherence_urea.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/singlet_states/decoherence_urea.m`
- Signature: `decoherence_urea()`
- Total lines: 52

## Purpose

A demonstration that the nitrogen singlet state in urea is not long-lived. The relaxation superoperator accounts for every di- polar coupling and every CSA tensor in the system. Calculation time: seconds

## Physical / mathematical content

- Long-lived singlet-state examples. The central concept is symmetry-protected or nearly symmetry-protected two-spin order that relaxes much more slowly than ordinary Zeeman magnetisation. Files here often analyse singlet-triplet subspaces, state conversion sequences, and relaxation leakage channels.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-14: Read the spin system (coordinates, chemical shifts, J-couplings and CSAs) from a vacuum DFT calculation; implemented by `[sys,inter]=g2spinach(gparse('../standard_systems/urea.log'), {{'H','1H'},{'N','15N'}},[30.0 166.0],[])`.
- Lines 16-17: Set magnet field to 1.0 Tesla; implemented by `sys.magnet=1.0`.
- Lines 19-20: Tighten up the tolerances; implemented by `sys.tols.rlx_integration=1e-5`.
- Lines 22-23: Set relaxation theory parameters; implemented by `inter.relaxation={'redfield'}`.
- Lines 28-29: Relaxation superoperator accuracy; implemented by `sys.tols.rlx_integration=1e-5`.
- Lines 32-33: Use complete basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 36-37: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 40-41: Build the relaxation superoperator; implemented by `R=relaxation(spin_system)`.
- Lines 43-44: Get the Lz and the singlet for 15N; implemented by `Lz=state(spin_system,'Lz','15N')`.
- Lines 47-48: See what remains of them under R; implemented by `disp(['Norm of R|Lz>: ' num2str(norm(R*Lz)/norm(Lz))])`.

### Key state/data transformations

- Lines 13-14: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../standard_systems/urea.log'), {{'H','1H'},{'N','15N'}},[30.0 166.0],[])`.
- Lines 17: computes `sys.magnet` using `sys.magnet=1.0`.
- Lines 20: computes `sys.tols.rlx_integration` using `sys.tols.rlx_integration=1e-5`.
- Lines 23: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 24: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 25: computes `inter.rlx_keep` using `inter.rlx_keep='labframe'`.
- Lines 26: computes `inter.tau_c` using `inter.tau_c={100e-12}`.
- Lines 30: computes `sys.tols.rlx_zero` using `sys.tols.rlx_zero=1e-5`.
- Lines 33: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 34: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 37: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 41: computes `R` using `R=relaxation(spin_system)`.
- Lines 44: computes `Lz` using `Lz=state(spin_system,'Lz','15N')`.
- Lines 45: computes `S` using `S=singlet(spin_system,1,4)`.

## Implementation structure

- A demonstration that the nitrogen singlet state in urea is not
- long-lived. The relaxation superoperator accounts for every di-
- polar coupling and every CSA tensor in the system.
- Calculation time: seconds
- Read the spin system (coordinates, chemical shifts,
- J-couplings and CSAs) from a vacuum DFT calculation
- Set magnet field to 1.0 Tesla
- Tighten up the tolerances
- Set relaxation theory parameters
- Relaxation superoperator accuracy
- Use complete basis set
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `relaxation()`, `state()`, `singlet()`, `num2str()`.
