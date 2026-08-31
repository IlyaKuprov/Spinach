# examples/quantum_tech/transmon_cavity_swap.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/transmon_cavity_swap.m`
- Signature: `transmon_cavity_swap()`
- Total lines: 65

## Purpose

Vacuum Rabi swap between a transmon and a microwave cavity mode, both represented by truncated bosonic Weyl algebras. This is the circuit-QED Jaynes-Cummings limit of Blais et al., Rev. Mod. Phys. 93, 025005 (2021). Calculation time: seconds

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- The effective hardware model is a weakly anharmonic oscillator. Duffing nonlinearity breaks equal level spacing and allows qubit-like addressability within a truncated bosonic ladder.
- The physics is Jaynes-Cummings-like cavity QED: a two-level or few-level matter degree of freedom exchanges excitations with a quantised harmonic mode through rotating terms such as a†σ_- + aσ_+.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Magnet field; implemented by `sys.magnet=0`.
- Lines 15-16: Particle specification; implemented by `sys.isotopes={'T3','C3'}`.
- Lines 18-19: Resonant transmon-cavity pair in the rotating frame; implemented by `inter.modes.frqs={0 0}`.
- Lines 24-25: Formalism and basis; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 28-29: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 32-33: Sequence parameters; implemented by `parameters.rho0=state(spin_system,{'BL2','BL1'},{1,2})`.
- Lines 37-38: Trajectory through the device context; implemented by `traj=device(spin_system,@traject,parameters,'cavity')`.
- Lines 40-41: Project out the transmon and cavity excitation populations; implemented by `coil_t=state(spin_system,{'BL2','E'},{1,2})`.
- Lines 46-47: Validate visible excitation exchange; implemented by `if (max(real(pop_c))<0.95)||(min(real(pop_t))>0.05)`.
- Lines 51-52: Validate population conservation in the active doublet; implemented by `if max(abs(real(pop_t+pop_c)-1))>1e-6`.
- Lines 56-57: Plot the swap dynamics; implemented by `time_axis=linspace(0,150,301)`.

### Control flow inferred from the code

- Line 47: conditional branch on `(max(real(pop_c))<0.95)||(min(real(pop_t))>0.05)`.
- Line 52: conditional branch on `max(abs(real(pop_t+pop_c)-1))>1e-6`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=0`.
- Lines 16: computes `sys.isotopes` using `sys.isotopes={'T3','C3'}`.
- Lines 19: computes `inter.modes.frqs` using `inter.modes.frqs={0 0}`.
- Lines 20: computes `inter.modes.anharms` using `inter.modes.anharms={-250e6 []}`.
- Lines 21: computes `inter.modes.exchange` using `inter.modes.exchange=cell(2,2)`.
- Lines 22: computes `inter.modes.exchange{1,2}` using `inter.modes.exchange{1,2}=20e6`.
- Lines 25: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 26: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 29: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 33: computes `parameters.rho0` using `parameters.rho0=state(spin_system,{'BL2','BL1'},{1,2})`.
- Lines 34: computes `parameters.sweep` using `parameters.sweep=2e9`.
- Lines 35: computes `parameters.npoints` using `parameters.npoints=301`.
- Lines 38: computes `traj` using `traj=device(spin_system,@traject,parameters,'cavity')`.
- Lines 41: computes `coil_t` using `coil_t=state(spin_system,{'BL2','E'},{1,2})`.
- Lines 42: computes `coil_c` using `coil_c=state(spin_system,{'E','BL2'},{1,2})`.
- Lines 43: computes `pop_t` using `pop_t=cellfun(@(rho)full(hdot(coil_t,rho)),traj)`.
- Lines 44: computes `pop_c` using `pop_c=cellfun(@(rho)full(hdot(coil_c,rho)),traj)`.
- Lines 57: computes `time_axis` using `time_axis=linspace(0,150,301)`.

## Implementation structure

- Vacuum Rabi swap between a transmon and a microwave cavity
- mode, both represented by truncated bosonic Weyl algebras.
- This is the circuit-QED Jaynes-Cummings limit of Blais et
- al., Rev. Mod. Phys. 93, 025005 (2021).
- Calculation time: seconds
- Magnet field
- Particle specification
- Resonant transmon-cavity pair in the rotating frame
- Formalism and basis
- Spinach housekeeping
- Sequence parameters
- Trajectory through the device context

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `device()`, `cellfun()`, `hdot()`, `kfigure()`, `kxlabel()`, `kylabel()`, `ktitle()`, `klegend()`.
