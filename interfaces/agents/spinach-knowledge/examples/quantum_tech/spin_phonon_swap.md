# examples/quantum_tech/spin_phonon_swap.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/spin_phonon_swap.m`
- Signature: `spin_phonon_swap()`
- Total lines: 64

## Purpose

Resonant excitation swap between an electron spin and a quantised phonon mode. The model is the spin-phonon Jaynes- Cummings limit used in mechanical spin-qubit proposals such as Rabl et al., Nature Physics 6, 602 (2010). Calculation time: seconds

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Magnet field; implemented by `sys.magnet=0`.
- Lines 15-16: Particle specification; implemented by `sys.isotopes={'E','V3'}`.
- Lines 18-19: Resonant phonon mode in the rotating frame; implemented by `inter.modes.frqs={[] 0}`.
- Lines 23-24: Formalism and basis; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 27-28: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 31-32: Sequence parameters; implemented by `parameters.rho0=state(spin_system,{'ZL2','BL1'},{1,2})`.
- Lines 36-37: Trajectory, 'cavity' is the set that keeps spin-mode exchange; implemented by `traj=device(spin_system,@traject,parameters,'cavity')`.
- Lines 39-40: Project out the spin and phonon excitation populations; implemented by `coil_s=state(spin_system,{'ZL2','E'},{1,2})`.
- Lines 45-46: Validate visible excitation exchange; implemented by `if (max(real(pop_v))<0.95)||(min(real(pop_s))>0.05)`.
- Lines 50-51: Validate population conservation in the active doublet; implemented by `if max(abs(real(pop_s+pop_v)-1))>1e-6`.
- Lines 55-56: Plot the swap dynamics; implemented by `time_axis=linspace(0,800,801)`.

### Control flow inferred from the code

- Line 46: conditional branch on `(max(real(pop_v))<0.95)||(min(real(pop_s))>0.05)`.
- Line 51: conditional branch on `max(abs(real(pop_s+pop_v)-1))>1e-6`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=0`.
- Lines 16: computes `sys.isotopes` using `sys.isotopes={'E','V3'}`.
- Lines 19: computes `inter.modes.frqs` using `inter.modes.frqs={[] 0}`.
- Lines 20: computes `inter.modes.exchange` using `inter.modes.exchange=cell(2,2)`.
- Lines 21: computes `inter.modes.exchange{1,2}` using `inter.modes.exchange{1,2}=4e6`.
- Lines 24: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 25: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 28: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 32: computes `parameters.rho0` using `parameters.rho0=state(spin_system,{'ZL2','BL1'},{1,2})`.
- Lines 33: computes `parameters.sweep` using `parameters.sweep=1e9`.
- Lines 34: computes `parameters.npoints` using `parameters.npoints=801`.
- Lines 37: computes `traj` using `traj=device(spin_system,@traject,parameters,'cavity')`.
- Lines 40: computes `coil_s` using `coil_s=state(spin_system,{'ZL2','E'},{1,2})`.
- Lines 41: computes `coil_v` using `coil_v=state(spin_system,{'ZL1','BL2'},{1,2})`.
- Lines 42: computes `pop_s` using `pop_s=cellfun(@(rho)full(hdot(coil_s,rho)),traj)`.
- Lines 43: computes `pop_v` using `pop_v=cellfun(@(rho)full(hdot(coil_v,rho)),traj)`.
- Lines 56: computes `time_axis` using `time_axis=linspace(0,800,801)`.

## Implementation structure

- Resonant excitation swap between an electron spin and a
- quantised phonon mode. The model is the spin-phonon Jaynes-
- Cummings limit used in mechanical spin-qubit proposals such
- as Rabl et al., Nature Physics 6, 602 (2010).
- Calculation time: seconds
- Magnet field
- Particle specification
- Resonant phonon mode in the rotating frame
- Formalism and basis
- Spinach housekeeping
- Sequence parameters
- Trajectory, 'cavity' is the set that keeps spin-mode exchange

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `device()`, `cellfun()`, `hdot()`, `kfigure()`, `ylim()`, `kxlabel()`, `kylabel()`, `ktitle()`, `klegend()`.
