# examples/singlet_states/decoherence_naphthalenetetrone.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/singlet_states/decoherence_naphthalenetetrone.m`
- Signature: `decoherence_naphthalenetetrone()`
- Total lines: 49

## Purpose

Long-lived spin states in the napthalenetetrone molecule. (4 protons, 256-dimensional Liouville space). The relaxation superoperator accounts for every dipolar coupling and every CSA tensor in the system. Calculation time: seconds

## Physical / mathematical content

- Long-lived singlet-state examples. The central concept is symmetry-protected or nearly symmetry-protected two-spin order that relaxes much more slowly than ordinary Zeeman magnetisation. Files here often analyse singlet-triplet subspaces, state conversion sequences, and relaxation leakage channels.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-16: Read the spin system (coordinates, chemical shifts, J-couplings and CSAs) from a vacuum DFT calculation; implemented by `[sys,inter]=g2spinach(gparse('../standard_systems/naphthalenetetrone.log'), {{'H','1H'}},31.8,[])`.
- Lines 17-18: Set magnet field to 1.0 Tesla; implemented by `sys.magnet=1.0`.
- Lines 20-21: Tighten up the tolerances; implemented by `sys.tols.rlx_integration=1e-5`.
- Lines 23-24: Set relaxation theory parameters; implemented by `inter.relaxation={'redfield'}`.
- Lines 29-30: Relaxation superoperator accuracy; implemented by `sys.tols.rlx_integration=1e-5`.
- Lines 33-34: Use complete basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 37-38: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 41-42: Build the relaxation superoperator; implemented by `R=relaxation(spin_system)`.
- Lines 44-45: List twenty smallest relaxation rates; implemented by `disp('Twenty smallest relaxation rates, Hz:')`.

### Key state/data transformations

- Lines 15-16: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../standard_systems/naphthalenetetrone.log'), {{'H','1H'}},31.8,[])`.
- Lines 18: computes `sys.magnet` using `sys.magnet=1.0`.
- Lines 21: computes `sys.tols.rlx_integration` using `sys.tols.rlx_integration=1e-5`.
- Lines 24: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 25: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 26: computes `inter.rlx_keep` using `inter.rlx_keep='labframe'`.
- Lines 27: computes `inter.tau_c` using `inter.tau_c={100e-12}`.
- Lines 31: computes `sys.tols.rlx_zero` using `sys.tols.rlx_zero=1e-5`.
- Lines 34: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 35: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 38: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 42: computes `R` using `R=relaxation(spin_system)`.

## Implementation structure

- Long-lived spin states in the napthalenetetrone molecule.
- (4 protons, 256-dimensional Liouville space). The relaxation
- superoperator accounts for every dipolar coupling and every
- CSA tensor in the system.
- Calculation time: seconds
- Read the spin system (coordinates, chemical shifts,
- J-couplings and CSAs) from a vacuum DFT calculation
- Set magnet field to 1.0 Tesla
- Tighten up the tolerances
- Set relaxation theory parameters
- Relaxation superoperator accuracy
- Use complete basis set

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `relaxation()`, `speye()`.
