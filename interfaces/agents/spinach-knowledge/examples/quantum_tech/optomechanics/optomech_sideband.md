# examples/quantum_tech/optomechanics/optomech_sideband.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/optomechanics/optomech_sideband.m`
- Signature: `optomech_sideband()`
- Total lines: 89

## Purpose

Optomechanical sideband transfer of a phonon Fock state into a driven cavity. A red-detuned coherent drive on the cavity acti- vates the beam-splitter part of the radiation pressure coupling, and the second Fock state of the mechanical oscillator is cohe- rently exchanged with the cavity field. Model and parameters from the propagation test set of QuantumPropagators.jl; all quantities are in the dimensionless units o

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Magnet field; implemented by `sys.magnet=0`.
- Lines 19-20: Cavity with five and phonon mode with eleven Fock levels; implemented by `sys.isotopes={'C5','V11'}`.
- Lines 22-23: Mode frequencies, cavity in the red-detuned drive rotating frame; implemented by `inter.modes.frqs={10/(2*pi) 10/(2*pi)}`.
- Lines 25-26: Radiation pressure coupling, cavity number times phonon coordinate; implemented by `inter.modes.longitudinal=cell(2,2)`.
- Lines 29-30: Formalism and basis; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 33-34: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 37-38: Hamiltonian from the declared interactions; implemented by `H=hamiltonian(assume(spin_system,'labframe'))`.
- Lines 40-41: Coherent drive on the cavity; implemented by `H=H+2*(operator(spin_system,'C',1)+operator(spin_system,'A',1))`.
- Lines 43-44: Number operators of the two modes; implemented by `Nc=operator(spin_system,'N',1)`.
- Lines 47-48: Empty cavity and second phonon Fock state; implemented by `rho=state(spin_system,{'BL1','BL3'},{1,2})`.
- Lines 50-51: Time grid of the source calculation; implemented by `nsteps=250; dt=0.2; time_axis=dt*(0:nsteps)`.
- Lines 53-54: Propagator over one time step; implemented by `P=propagator(spin_system,H,dt)`.
- Lines 56-57: Preallocate occupation trajectories; implemented by `n_cav=zeros(nsteps+1,1); n_mech=zeros(nsteps+1,1)`.
- Lines 59-60: Propagate and record mode occupations; implemented by `for n=1:(nsteps+1)`.
- Lines 66-67: Validate the trace preservation; implemented by `if abs(trace(rho)-1)>1e-6`.
- Lines 71-72: Validate the initial phonon occupation; implemented by `if abs(n_mech(1)-2)>1e-12`.
- Lines 76-77: Validate the sideband transfer into the cavity; implemented by `if max(n_cav)<1.5`.
- Lines 81-82: Plot the mode occupation dynamics; implemented by `kfigure(); plot(time_axis,[n_cav n_mech],'LineWidth',1.5)`.

### Control flow inferred from the code

- Line 60: `for` loop over `n=1:(nsteps+1)`.
- Line 67: conditional branch on `abs(trace(rho)-1)>1e-6`.
- Line 72: conditional branch on `abs(n_mech(1)-2)>1e-12`.
- Line 77: conditional branch on `max(n_cav)<1.5`.

### Key state/data transformations

- Lines 17: computes `sys.magnet` using `sys.magnet=0`.
- Lines 20: computes `sys.isotopes` using `sys.isotopes={'C5','V11'}`.
- Lines 23: computes `inter.modes.frqs` using `inter.modes.frqs={10/(2*pi) 10/(2*pi)}`.
- Lines 26: computes `inter.modes.longitudinal` using `inter.modes.longitudinal=cell(2,2)`.
- Lines 27: computes `inter.modes.longitudinal{1,2}` using `inter.modes.longitudinal{1,2}=-sqrt(2)/(2*pi)`.
- Lines 30: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 31: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 34: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 38: computes `H` using `H=hamiltonian(assume(spin_system,'labframe'))`.
- Lines 44: computes `Nc` using `Nc=operator(spin_system,'N',1)`.
- Lines 45: computes `Nm` using `Nm=operator(spin_system,'N',2)`.
- Lines 48: computes `rho` using `rho=state(spin_system,{'BL1','BL3'},{1,2})`.
- Lines 51: computes `nsteps` using `nsteps=250; dt=0.2; time_axis=dt*(0:nsteps)`.
- Lines 54: computes `P` using `P=propagator(spin_system,H,dt)`.
- Lines 57: computes `n_cav` using `n_cav=zeros(nsteps+1,1); n_mech=zeros(nsteps+1,1)`.
- Lines 61: computes `n_cav(n)` using `n_cav(n)=real(trace(Nc*rho))`.
- Lines 62: computes `n_mech(n)` using `n_mech(n)=real(trace(Nm*rho))`.

## Implementation structure

- Optomechanical sideband transfer of a phonon Fock state into a
- driven cavity. A red-detuned coherent drive on the cavity acti-
- vates the beam-splitter part of the radiation pressure coupling,
- and the second Fock state of the mechanical oscillator is cohe-
- rently exchanged with the cavity field. Model and parameters
- from the propagation test set of QuantumPropagators.jl; all
- quantities are in the dimensionless units of the source, with
- input values divided by 2*pi to cancel the Hz convention.
- Calculation time: seconds
- Magnet field
- Cavity with five and phonon mode with eleven Fock levels
- Mode frequencies, cavity in the red-detuned drive rotating frame

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `hamiltonian()`, `assume()`, `operator()`, `state()`, `propagator()`, `n_cav()`, `n_mech()`, `kfigure()`, `kxlabel()`, `kylabel()`, `ktitle()`, `klegend()`.
