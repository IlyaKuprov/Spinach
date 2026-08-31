# examples/quantum_tech/spin_cavity_vacuum_rabi.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/spin_cavity_vacuum_rabi.m`
- Signature: `spin_cavity_vacuum_rabi()`
- Total lines: 65

## Purpose

Vacuum Rabi oscillation between an electron spin and a micro- wave cavity mode in the Jaynes-Cummings approximation. This is the one-spin limit of the spin-ensemble cavity experiments of Schuster et al. and Kubo et al., Phys. Rev. Lett. 105, 140501 and 140502 (2010). Calculation time: seconds

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- The physics is Jaynes-Cummings-like cavity QED: a two-level or few-level matter degree of freedom exchanges excitations with a quantised harmonic mode through rotating terms such as a†σ_- + aσ_+.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Magnet field; implemented by `sys.magnet=0`.
- Lines 16-17: Particle specification; implemented by `sys.isotopes={'E','C3'}`.
- Lines 19-20: Resonant cavity in the rotating frame; implemented by `inter.modes.frqs={[] 0}`.
- Lines 24-25: Formalism and basis; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 28-29: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 32-33: Sequence parameters; implemented by `parameters.rho0=state(spin_system,{'ZL2','BL1'},{1,2})`.
- Lines 37-38: Trajectory through the device context; implemented by `traj=device(spin_system,@traject,parameters,'cavity')`.
- Lines 40-41: Project out the spin and cavity excitation populations; implemented by `coil_s=state(spin_system,{'ZL2','E'},{1,2})`.
- Lines 46-47: Validate visible excitation exchange; implemented by `if (max(real(pop_c))<0.95)||(min(real(pop_s))>0.05)`.
- Lines 51-52: Validate population conservation in the active doublet; implemented by `if max(abs(real(pop_s+pop_c)-1))>1e-6`.
- Lines 56-57: Plot the vacuum Rabi dynamics; implemented by `time_axis=linspace(0,500,501)`.

### Control flow inferred from the code

- Line 47: conditional branch on `(max(real(pop_c))<0.95)||(min(real(pop_s))>0.05)`.
- Line 52: conditional branch on `max(abs(real(pop_s+pop_c)-1))>1e-6`.

### Key state/data transformations

- Lines 14: computes `sys.magnet` using `sys.magnet=0`.
- Lines 17: computes `sys.isotopes` using `sys.isotopes={'E','C3'}`.
- Lines 20: computes `inter.modes.frqs` using `inter.modes.frqs={[] 0}`.
- Lines 21: computes `inter.modes.exchange` using `inter.modes.exchange=cell(2,2)`.
- Lines 22: computes `inter.modes.exchange{1,2}` using `inter.modes.exchange{1,2}=8e6`.
- Lines 25: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 26: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 29: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 33: computes `parameters.rho0` using `parameters.rho0=state(spin_system,{'ZL2','BL1'},{1,2})`.
- Lines 34: computes `parameters.sweep` using `parameters.sweep=1e9`.
- Lines 35: computes `parameters.npoints` using `parameters.npoints=501`.
- Lines 38: computes `traj` using `traj=device(spin_system,@traject,parameters,'cavity')`.
- Lines 41: computes `coil_s` using `coil_s=state(spin_system,{'ZL2','E'},{1,2})`.
- Lines 42: computes `coil_c` using `coil_c=state(spin_system,{'ZL1','BL2'},{1,2})`.
- Lines 43: computes `pop_s` using `pop_s=cellfun(@(rho)full(hdot(coil_s,rho)),traj)`.
- Lines 44: computes `pop_c` using `pop_c=cellfun(@(rho)full(hdot(coil_c,rho)),traj)`.
- Lines 57: computes `time_axis` using `time_axis=linspace(0,500,501)`.

## Implementation structure

- Vacuum Rabi oscillation between an electron spin and a micro-
- wave cavity mode in the Jaynes-Cummings approximation. This
- is the one-spin limit of the spin-ensemble cavity experiments
- of Schuster et al. and Kubo et al., Phys. Rev. Lett. 105,
- 140501 and 140502 (2010).
- Calculation time: seconds
- Magnet field
- Particle specification
- Resonant cavity in the rotating frame
- Formalism and basis
- Spinach housekeeping
- Sequence parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `device()`, `cellfun()`, `hdot()`, `kfigure()`, `ylim()`, `kxlabel()`, `kylabel()`, `ktitle()`, `klegend()`.
