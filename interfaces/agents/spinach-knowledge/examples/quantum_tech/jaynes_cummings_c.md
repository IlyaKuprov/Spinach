# examples/quantum_tech/jaynes_cummings_c.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/jaynes_cummings_c.m`
- Signature: `jaynes_cummings_c()`
- Total lines: 74

## Purpose

An exchange-coupled two-electron system with the electrons having independent Jaynes-Cummings couplings to the same mode of an electromagnetic cavity. A time-domain simulati- on starting with transverse spin magnetisation and empty cavity mode. Detected on the Lx operator of the spin and magnetic field operator of the cavity mode. Calculation time: seconds

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- The physics is Jaynes-Cummings-like cavity QED: a two-level or few-level matter degree of freedom exchanges excitations with a quantised harmonic mode through rotating terms such as a†σ_- + aσ_+.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: Magnet field; implemented by `sys.magnet=0.33`.
- Lines 17-18: System; implemented by `sys.isotopes={'E','E','C5'}`.
- Lines 20-21: Exchange coupling between the electrons; implemented by `inter.coupling.scalar=cell(3,3)`.
- Lines 24-25: Cavity resonant with the electrons; implemented by `e_frq=-sys.magnet*spin('E')/(2*pi)`.
- Lines 31-32: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 35-36: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 39-40: Sequence parameters; implemented by `parameters.spins={'E'}`.
- Lines 47-48: Trajectory through the device context; implemented by `traj=device(spin_system,@traject,parameters,'cavity')`.
- Lines 50-52: Detection state, spins; implemented by `S=state(spin_system,{'Lx'},{1})+ state(spin_system,{'Lx'},{2})`.
- Lines 54-56: Detection state, cavity mode; implemented by `B=(state(spin_system,{'C'},{3})- state(spin_system,{'A'},{3}))/2i`.
- Lines 58-59: Project out the observables; implemented by `traj_s=S'*traj; traj_c=B'*traj`.
- Lines 61-62: Plot the observables; implemented by `time_axis=linspace(0,2.5,251)`.

### Key state/data transformations

- Lines 15: computes `sys.magnet` using `sys.magnet=0.33`.
- Lines 18: computes `sys.isotopes` using `sys.isotopes={'E','E','C5'}`.
- Lines 21: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(3,3)`.
- Lines 22: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=5e6`.
- Lines 25: computes `e_frq` using `e_frq=-sys.magnet*spin('E')/(2*pi)`.
- Lines 26: computes `inter.modes.frqs` using `inter.modes.frqs={[] [] e_frq}`.
- Lines 27: computes `inter.modes.exchange` using `inter.modes.exchange=cell(3,3)`.
- Lines 28: computes `inter.modes.exchange{1,3}` using `inter.modes.exchange{1,3}=2.828e6`.
- Lines 29: computes `inter.modes.exchange{2,3}` using `inter.modes.exchange{2,3}=2.728e6`.
- Lines 32: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 33: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 36: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 40: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 41: computes `parameters.offset` using `parameters.offset=5e6`.
- Lines 42-43: computes `parameters.rho0` using `parameters.rho0=state(spin_system,{'Lx','BL1'},{1,3})+ state(spin_system,{'Lx','BL1'},{2,3})`.
- Lines 44: computes `parameters.sweep` using `parameters.sweep=1e8`.
- Lines 45: computes `parameters.npoints` using `parameters.npoints=251`.
- Lines 48: computes `traj` using `traj=device(spin_system,@traject,parameters,'cavity')`.

## Implementation structure

- An exchange-coupled two-electron system with the electrons
- having independent Jaynes-Cummings couplings to the same
- mode of an electromagnetic cavity. A time-domain simulati-
- on starting with transverse spin magnetisation and empty
- cavity mode. Detected on the Lx operator of the spin and
- magnetic field operator of the cavity mode.
- Calculation time: seconds
- Magnet field
- System
- Exchange coupling between the electrons
- Cavity resonant with the electrons
- Basis set

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `spin()`, `create()`, `basis()`, `state()`, `device()`, `kfigure()`, `scale_figure()`, `subplot()`, `kxlabel()`, `kylabel()`, `ktitle()`.
