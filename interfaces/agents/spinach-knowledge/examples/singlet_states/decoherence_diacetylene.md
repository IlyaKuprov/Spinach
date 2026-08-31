# examples/singlet_states/decoherence_diacetylene.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/singlet_states/decoherence_diacetylene.m`
- Signature: `decoherence_diacetylene()`
- Total lines: 67

## Purpose

Long-lived spin states in the diacetylene molecule (2 protons, 4 carbons, 4096-dimensional Liouville space). The relaxation superoperator accounts for every dipolar coupling and every CSA tensor in the system. Calculation time: seconds

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

- Lines 13-16: Read the spin system (coordinates, chemical shifts, J-couplings and CSAs) from a vacuum DFT calculation; implemented by `[sys,inter]=g2spinach(gparse('../standard_systems/diacetylene.log'), {{'H','1H'},{'C','13C'}},[31.8 182.4],[])`.
- Lines 18-19: Set magnet field to 1.0 Tesla; implemented by `sys.magnet=14.1`.
- Lines 21-22: Tighten up the tolerances; implemented by `sys.tols.rlx_integration=1e-5`.
- Lines 24-25: Set relaxation theory parameters; implemented by `inter.relaxation={'redfield'}`.
- Lines 30-31: Relaxation superoperator accuracy; implemented by `sys.tols.rlx_integration=1e-5`.
- Lines 34-35: Use complete basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 38-39: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 42-43: Build the relaxation superoperator; implemented by `R=relaxation(spin_system)`.
- Lines 46-48: List twenty smallest magnitude eigenvalues of the relaxation superoperator (diagonal preconditioning is used); implemented by `disp('20 smallest eigenvalues of the relaxation superoperator:')`.
- Lines 51-53: Compute the self-relaxation rate of the singlet state between the two centre carbons (spins 1 and 2 in this case); implemented by `S=singlet(spin_system,1,2); S=S/norm(S)`.
- Lines 57-58: Find the eigenvectors corresponding to the slowly relaxing states; implemented by `[v,~]=eigs(R-speye(size(R)),2,'SM')`.
- Lines 60-61: Get the spherical tensor composition of the slowly relaxing states; implemented by `disp('Slowly relaxing state 1:')`.

### Key state/data transformations

- Lines 15-16: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../standard_systems/diacetylene.log'), {{'H','1H'},{'C','13C'}},[31.8 182.4],[])`.
- Lines 19: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 22: computes `sys.tols.rlx_integration` using `sys.tols.rlx_integration=1e-5`.
- Lines 25: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 26: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 27: computes `inter.rlx_keep` using `inter.rlx_keep='labframe'`.
- Lines 28: computes `inter.tau_c` using `inter.tau_c={100e-12}`.
- Lines 32: computes `sys.tols.rlx_zero` using `sys.tols.rlx_zero=1e-5`.
- Lines 35: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 36: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 39: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 43: computes `R` using `R=relaxation(spin_system)`.
- Lines 53: computes `S` using `S=singlet(spin_system,1,2); S=S/norm(S)`.
- Lines 58: computes `[v,~]` using `[v,~]=eigs(R-speye(size(R)),2,'SM')`.

## Implementation structure

- Long-lived spin states in the diacetylene molecule (2 protons,
- 4 carbons, 4096-dimensional Liouville space). The relaxation
- superoperator accounts for every dipolar coupling and every CSA
- tensor in the system.
- Calculation time: seconds
- Read the spin system (coordinates, chemical shifts,
- J-couplings and CSAs) from a vacuum DFT calculation
- Set magnet field to 1.0 Tesla
- Tighten up the tolerances
- Set relaxation theory parameters
- Relaxation superoperator accuracy
- Use complete basis set

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `relaxation()`, `format()`, `speye()`, `singlet()`, `stateinfo()`.
