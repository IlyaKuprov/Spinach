# examples/quantum_tech/spin_phonon_dephasing.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/spin_phonon_dephasing.m`
- Signature: `spin_phonon_dephasing()`
- Total lines: 72

## Purpose

Longitudinal spin-phonon coupling producing spin coherence modulation and spin-conditioned displacement of a quantised vibrational mode. This is a minimal Weyl-algebra version of strain-modulated spin Hamiltonians used for NV-centre and molecular spin-phonon dynamics. Both observables come from a single trajectory: the drift commutes with the spin projec- tion operator, so the coherence and the population sectors of 

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Magnet field; implemented by `sys.magnet=0`.
- Lines 20-21: Particle specification; implemented by `sys.isotopes={'E','V7'}`.
- Lines 23-24: Phonon mode with a longitudinal coupling; implemented by `inter.modes.frqs={[] 20e6}`.
- Lines 28-29: Formalism and basis; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 32-33: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 36-37: Sequence parameters; implemented by `parameters.sweep=5e8`.
- Lines 42-43: Trajectory through the device context; implemented by `traj=device(spin_system,@traject,parameters,'spin-phonon')`.
- Lines 45-46: Project out the spin coherence and the conditional displacement; implemented by `coil_s=state(spin_system,'Lx',1)`.
- Lines 51-52: Validate visible oscillator displacement; implemented by `if max(abs(real(traj_q)))<0.2`.
- Lines 56-57: Validate visible spin coherence modulation; implemented by `if (max(real(traj_s))-min(real(traj_s)))<0.02`.
- Lines 61-62: Plot the coupled dynamics; implemented by `time_axis=linspace(0,1.0,501)`.

### Control flow inferred from the code

- Line 52: conditional branch on `max(abs(real(traj_q)))<0.2`.
- Line 57: conditional branch on `(max(real(traj_s))-min(real(traj_s)))<0.02`.

### Key state/data transformations

- Lines 18: computes `sys.magnet` using `sys.magnet=0`.
- Lines 21: computes `sys.isotopes` using `sys.isotopes={'E','V7'}`.
- Lines 24: computes `inter.modes.frqs` using `inter.modes.frqs={[] 20e6}`.
- Lines 25: computes `inter.modes.longitudinal` using `inter.modes.longitudinal=cell(2,2)`.
- Lines 26: computes `inter.modes.longitudinal{1,2}` using `inter.modes.longitudinal{1,2}=4e6*sqrt(2)`.
- Lines 29: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 30: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 33: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 37: computes `parameters.sweep` using `parameters.sweep=5e8`.
- Lines 38: computes `parameters.npoints` using `parameters.npoints=501`.
- Lines 39-40: computes `parameters.rho0` using `parameters.rho0=state(spin_system,{'Lx','BL1'},{1,2})+ state(spin_system,{'ZL2','BL1'},{1,2})`.
- Lines 43: computes `traj` using `traj=device(spin_system,@traject,parameters,'spin-phonon')`.
- Lines 46: computes `coil_s` using `coil_s=state(spin_system,'Lx',1)`.
- Lines 47: computes `coil_q` using `coil_q=state(spin_system,'C',2)+state(spin_system,'A',2)`.
- Lines 48: computes `traj_s` using `traj_s=cellfun(@(rho)full(hdot(coil_s,rho)),traj)`.
- Lines 49: computes `traj_q` using `traj_q=cellfun(@(rho)full(hdot(coil_q,rho)),traj)`.
- Lines 62: computes `time_axis` using `time_axis=linspace(0,1.0,501)`.

## Implementation structure

- Longitudinal spin-phonon coupling producing spin coherence
- modulation and spin-conditioned displacement of a quantised
- vibrational mode. This is a minimal Weyl-algebra version of
- strain-modulated spin Hamiltonians used for NV-centre and
- molecular spin-phonon dynamics. Both observables come from a
- single trajectory: the drift commutes with the spin projec-
- tion operator, so the coherence and the population sectors
- of the initial condition evolve and are detected without
- mixing.
- Calculation time: seconds
- Magnet field
- Particle specification

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `device()`, `cellfun()`, `hdot()`, `kfigure()`, `scale_figure()`, `subplot()`, `kxlabel()`, `kylabel()`, `ktitle()`.
