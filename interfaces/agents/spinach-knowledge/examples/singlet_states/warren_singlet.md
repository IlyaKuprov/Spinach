# examples/singlet_states/warren_singlet.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/singlet_states/warren_singlet.m`
- Signature: `warren_singlet()`
- Total lines: 45

## Purpose

A demonstration that long-lived states exist that are immune not only to dipolar and CSA, but also to quadrupolar relaxati- on in certain circumstances. Full Redfield superoperator for dipolar and quadrupolar relaxation in liquid state is compu- ted and diagonalized. Calculation time: seconds

## Physical / mathematical content

- Long-lived singlet-state examples. The central concept is symmetry-protected or nearly symmetry-protected two-spin order that relaxes much more slowly than ordinary Zeeman magnetisation. Files here often analyse singlet-triplet subspaces, state conversion sequences, and relaxation leakage channels.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: System specification; implemented by `sys.magnet=14.1`.
- Lines 22-23: Relaxation theory parameters; implemented by `inter.relaxation={'redfield'}`.
- Lines 28-29: Relaxation superoperator accuracy; implemented by `sys.tols.rlx_integration=1e-5`.
- Lines 32-33: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 36-37: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 40-41: Relaxation superoperator; implemented by `R=relaxation(spin_system)`.

### Key state/data transformations

- Lines 15: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 16: computes `sys.isotopes` using `sys.isotopes={'14N','14N'}`.
- Lines 17: computes `inter.coordinates` using `inter.coordinates={[0.0 0.0 0.0]`.
- Lines 19: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=eeqq2nqi(1.25e6,0.25,1,[0 0 0])`.
- Lines 20: computes `inter.coupling.matrix{2,2}` using `inter.coupling.matrix{2,2}=eeqq2nqi(1.25e6,0.25,1,[0 0 0])`.
- Lines 23: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 24: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 25: computes `inter.rlx_keep` using `inter.rlx_keep='labframe'`.
- Lines 26: computes `inter.tau_c` using `inter.tau_c={5e-9}`.
- Lines 29: computes `sys.tols.rlx_integration` using `sys.tols.rlx_integration=1e-5`.
- Lines 30: computes `sys.tols.rlx_zero` using `sys.tols.rlx_zero=1e-5`.
- Lines 33: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 34: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 37: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 41: computes `R` using `R=relaxation(spin_system)`.

## Implementation structure

- A demonstration that long-lived states exist that are immune
- not only to dipolar and CSA, but also to quadrupolar relaxati-
- on in certain circumstances. Full Redfield superoperator for
- dipolar and quadrupolar relaxation in liquid state is compu-
- ted and diagonalized.
- Calculation time: seconds
- System specification
- Relaxation theory parameters
- Relaxation superoperator accuracy
- Basis set
- Spinach housekeeping
- Relaxation superoperator

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `relaxation()`.
