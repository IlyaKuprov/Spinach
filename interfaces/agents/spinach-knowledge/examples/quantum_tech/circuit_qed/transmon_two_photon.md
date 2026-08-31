# examples/quantum_tech/circuit_qed/transmon_two_photon.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/circuit_qed/transmon_two_photon.m`
- Signature: `transmon_two_photon()`
- Total lines: 95

## Purpose

Two-photon transition in a four-level Duffing transmon. The drive carrier is placed halfway between the 0-1 and 1-2 tran- sition frequencies, where neither single-photon transition is resonant, and GRAPE finds a pulse that moves the population from the ground state into the second excited state through a virtual intermediate state. Model and parameters from Example 2 of the GRAPE_SCQ package. Calculation time: minute

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 15-16: Magnet field; implemented by `sys.magnet=0`.
- Lines 18-19: Particle specification; implemented by `sys.isotopes={'T4'}`.
- Lines 21-22: Transmon detuning from the carrier and anharmonicity; implemented by `inter.modes.frqs={100e6}`.
- Lines 25-26: Formalism and basis; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 29-30: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 33-34: Drift Hamiltonian from the declared interactions; implemented by `H=hamiltonian(assume(spin_system,'cavity'))`.
- Lines 36-37: Quadrature control operators of the transmon mode; implemented by `Cr=operator(spin_system,'C',1); An=operator(spin_system,'A',1)`.
- Lines 40-41: Build source and destination states; implemented by `rho_init=state(spin_system,'BL1',1)`.
- Lines 46-47: Define control parameters; implemented by `control.isotopes={'T4'}`.
- Lines 60-61: Plots during optimisation; implemented by `control.plotting={'xy_controls'}`.
- Lines 63-64: Spinach housekeeping; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 66-67: Smooth deterministic initial guess; implemented by `guess=[0.1*sin(pi*(1:200)/200); 0.3*cos(pi*(1:200)/200)]`.
- Lines 69-70: Run the optimisation, get normalised pulse; implemented by `pulse=fmaxnewton(spin_system,@grape_xy,guess)`.
- Lines 72-73: Propagate the optimal pulse and record populations; implemented by `dts=spin_system.control.pulse_dt; nslices=numel(dts)`.
- Lines 81-82: Validate the two-photon transfer; implemented by `if pops(3,end)<0.99`.
- Lines 86-87: Plot the level populations under the optimal pulse; implemented by `kfigure(); plot(1e9*[0 cumsum(dts)],pops','LineWidth',1.5)`.

### Control flow inferred from the code

- Line 75: `for` loop over `n=1:nslices`.
- Line 82: conditional branch on `pops(3,end)<0.99`.

### Key state/data transformations

- Lines 16: computes `sys.magnet` using `sys.magnet=0`.
- Lines 19: computes `sys.isotopes` using `sys.isotopes={'T4'}`.
- Lines 22: computes `inter.modes.frqs` using `inter.modes.frqs={100e6}`.
- Lines 23: computes `inter.modes.anharms` using `inter.modes.anharms={-200e6}`.
- Lines 26: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 27: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 30: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 34: computes `H` using `H=hamiltonian(assume(spin_system,'cavity'))`.
- Lines 37: computes `Cr` using `Cr=operator(spin_system,'C',1); An=operator(spin_system,'A',1)`.
- Lines 38: computes `Cx` using `Cx=(Cr+An)/2; Cy=1i*(Cr-An)/2`.
- Lines 41: computes `rho_init` using `rho_init=state(spin_system,'BL1',1)`.
- Lines 42: computes `rho_targ` using `rho_targ=state(spin_system,'BL3',1)`.
- Lines 47: computes `control.isotopes` using `control.isotopes={'T4'}`.
- Lines 48: computes `control.channels` using `control.channels=[1;1]`.
- Lines 49: computes `control.drifts` using `control.drifts={{H}}`.
- Lines 50: computes `control.operators` using `control.operators={Cx,Cy}`.
- Lines 51: computes `control.rho_init` using `control.rho_init={rho_init}`.
- Lines 52: computes `control.rho_targ` using `control.rho_targ={rho_targ}`.

## Implementation structure

- Two-photon transition in a four-level Duffing transmon. The
- drive carrier is placed halfway between the 0-1 and 1-2 tran-
- sition frequencies, where neither single-photon transition is
- resonant, and GRAPE finds a pulse that moves the population
- from the ground state into the second excited state through a
- virtual intermediate state. Model and parameters from Example
- 2 of the GRAPE_SCQ package.
- Calculation time: minutes
- Magnet field
- Particle specification
- Transmon detuning from the carrier and anharmonicity
- Formalism and basis

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `hamiltonian()`, `assume()`, `operator()`, `state()`, `optimcon()`, `fmaxnewton()`, `pops()`, `pulse()`, `propagator()`, `dts()`, `kfigure()`, `cumsum()`, `kxlabel()`, `kylabel()`.
