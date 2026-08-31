# examples/quantum_tech/jaynes_cummings_b.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/jaynes_cummings_b.m`
- Signature: `jaynes_cummings_b()`
- Total lines: 79

## Purpose

Jaynes-Cummings coupling between a spin and an electromagnetic cavity mode with five population numbers included. A time-dom- ain simulation starting with transverse spin magnetisation and empty cavity mode. Detected on the Lx operator of the spin and magnetic field operator of the cavity mode. Calculation time: seconds

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- The physics is Jaynes-Cummings-like cavity QED: a two-level or few-level matter degree of freedom exchanges excitations with a quantised harmonic mode through rotating terms such as a†σ_- + aσ_+.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Magnet field; implemented by `sys.magnet=0.33`.
- Lines 16-17: System; implemented by `sys.isotopes={'E','C5'}`.
- Lines 19-20: Cavity resonant with the electron; implemented by `e_frq=-sys.magnet*spin('E')/(2*pi)`.
- Lines 25-26: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 29-30: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 33-34: Sequence parameters; implemented by `parameters.spins={'E'}`.
- Lines 41-42: Trajectory through the device context; implemented by `traj=device(spin_system,@traject,parameters,'cavity')`.
- Lines 44-45: Detection state, spin; implemented by `S=state(spin_system,{'Lx'},{1})`.
- Lines 47-49: Detection state, cavity mode; implemented by `B=(state(spin_system,{'C'},{2})- state(spin_system,{'A'},{2}))/2i`.
- Lines 51-52: Project out the observables; implemented by `traj_s=S'*traj; traj_c=B'*traj`.
- Lines 54-55: Plot the observables; implemented by `time_axis=linspace(0,2.5,251)`.
- Lines 66-67: Cavity energy level population evolution; implemented by `L1=state(spin_system,{'E','BL1'},{1,2})`.
- Lines 72-73: Plot the level populations; implemented by `kfigure(); plot(time_axis,real([pop_1; pop_2; pop_3]))`.

### Key state/data transformations

- Lines 14: computes `sys.magnet` using `sys.magnet=0.33`.
- Lines 17: computes `sys.isotopes` using `sys.isotopes={'E','C5'}`.
- Lines 20: computes `e_frq` using `e_frq=-sys.magnet*spin('E')/(2*pi)`.
- Lines 21: computes `inter.modes.frqs` using `inter.modes.frqs={[] e_frq}`.
- Lines 22: computes `inter.modes.exchange` using `inter.modes.exchange=cell(2,2)`.
- Lines 23: computes `inter.modes.exchange{1,2}` using `inter.modes.exchange{1,2}=2.828e6`.
- Lines 26: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 27: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 30: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 34: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 35: computes `parameters.offset` using `parameters.offset=5e6`.
- Lines 36-37: computes `parameters.rho0` using `parameters.rho0=state(spin_system,{'Lx','E' },{1,2})+ state(spin_system,{'E' ,'BL1'},{1,2})/2`.
- Lines 38: computes `parameters.sweep` using `parameters.sweep=1e8`.
- Lines 39: computes `parameters.npoints` using `parameters.npoints=251`.
- Lines 42: computes `traj` using `traj=device(spin_system,@traject,parameters,'cavity')`.
- Lines 45: computes `S` using `S=state(spin_system,{'Lx'},{1})`.
- Lines 48-49: computes `B` using `B=(state(spin_system,{'C'},{2})- state(spin_system,{'A'},{2}))/2i`.
- Lines 52: computes `traj_s` using `traj_s=S'*traj; traj_c=B'*traj`.

## Implementation structure

- Jaynes-Cummings coupling between a spin and an electromagnetic
- cavity mode with five population numbers included. A time-dom-
- ain simulation starting with transverse spin magnetisation and
- empty cavity mode. Detected on the Lx operator of the spin and
- magnetic field operator of the cavity mode.
- Calculation time: seconds
- Magnet field
- System
- Cavity resonant with the electron
- Basis set
- Spinach housekeeping
- Sequence parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `spin()`, `create()`, `basis()`, `state()`, `device()`, `kfigure()`, `scale_figure()`, `subplot()`, `kxlabel()`, `kylabel()`, `ktitle()`, `klegend()`.
