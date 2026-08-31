# examples/singlet_states/m2s_example.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/singlet_states/m2s_example.m`
- Signature: `m2s_example()`
- Total lines: 45

## Purpose

An example of the M2S sequence for a two-spin system. Calculation time: seconds

## Physical / mathematical content

- Long-lived singlet-state examples. The central concept is symmetry-protected or nearly symmetry-protected two-spin order that relaxes much more slowly than ordinary Zeeman magnetisation. Files here often analyse singlet-triplet subspaces, state conversion sequences, and relaxation leakage channels.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Spin system and interactions; implemented by `sys.magnet=9.4`.
- Lines 17-18: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 21-22: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 25-26: Hamiltonian; implemented by `H=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 28-29: Pulse operators; implemented by `Cx=operator(spin_system,'Lx','13C')`.
- Lines 32-33: Start with longitudinal magnetisation; implemented by `rho0=state(spin_system,'Lz','all')`.
- Lines 35-36: Detect singlet state; implemented by `coil=singlet(spin_system,1,2)`.
- Lines 38-39: Call the M2S sequence; implemented by `rho=m2s(spin_system,H,Cx,Cy,rho0,55,6.0)`.
- Lines 41-42: Display the singlet population; implemented by `disp(['Singlet population: ' num2str(coil'*rho)])`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=9.4`.
- Lines 12: computes `sys.isotopes` using `sys.isotopes={'13C','13C'}`.
- Lines 13: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.03,-0.03}`.
- Lines 14: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(2)`.
- Lines 15: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=55`.
- Lines 18: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 19: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 22: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 26: computes `H` using `H=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 29: computes `Cx` using `Cx=operator(spin_system,'Lx','13C')`.
- Lines 30: computes `Cy` using `Cy=operator(spin_system,'Ly','13C')`.
- Lines 33: computes `rho0` using `rho0=state(spin_system,'Lz','all')`.
- Lines 36: computes `coil` using `coil=singlet(spin_system,1,2)`.
- Lines 39: computes `rho` using `rho=m2s(spin_system,H,Cx,Cy,rho0,55,6.0)`.

## Implementation structure

- An example of the M2S sequence for a two-spin system.
- Calculation time: seconds
- Spin system and interactions
- Basis set
- Spinach housekeeping
- Hamiltonian
- Pulse operators
- Start with longitudinal magnetisation
- Detect singlet state
- Call the M2S sequence
- Display the singlet population

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `hamiltonian()`, `assume()`, `operator()`, `state()`, `singlet()`, `m2s()`, `num2str()`.
