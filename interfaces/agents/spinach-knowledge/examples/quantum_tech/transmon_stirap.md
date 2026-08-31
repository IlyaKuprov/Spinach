# examples/quantum_tech/transmon_stirap.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/transmon_stirap.m`
- Signature: `transmon_stirap()`
- Total lines: 85

## Purpose

Basic single transmon system with Duffing model in- teractions, parameters and model from: Ensemble GRAPE optimisation with a distribution of control powers and a penalty on excess power. Calculation time: minutes

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.
- The effective hardware model is a weakly anharmonic oscillator. Duffing nonlinearity breaks equal level spacing and allows qubit-like addressability within a truncated bosonic ladder.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 15-16: Magnet field; implemented by `sys.magnet=0`.
- Lines 18-19: Particle specification; implemented by `sys.isotopes={'T3'}`.
- Lines 21-22: Rotating frame ladder detunings; implemented by `inter.modes.frqs={10e6}`.
- Lines 25-26: Formalism and basis; implemented by `bas.formalism='zeeman-liouv'`.
- Lines 29-30: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 33-34: Drift Hamiltonian from the declared interactions; implemented by `H=hamiltonian(assume(spin_system,'cavity'))`.
- Lines 36-37: Pulse power ensemble; implemented by `pwr_levels=2*pi*[30e6 35e6 40e6 45e6 50e6]`.
- Lines 39-40: Gaussian STIRAP pulse pair, normalised to the top power level; implemented by `t=1e-9*linspace(-150,150,300)`.
- Lines 47-48: Build control operators; implemented by `H01=sparse([0 1 0; 1 0 0; 0 0 0]/2)`.
- Lines 51-52: Build source and destination states; implemented by `rho_init=state(spin_system,'BL1',1)`.
- Lines 57-58: Define control parameters; implemented by `control.isotopes={'T3'}`.
- Lines 72-73: Plots during optimisation; implemented by `control.plotting={'xy_controls','spectrogram','robustness'}`.
- Lines 75-76: Spinach housekeeping; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 78-79: Gaussian pulse pair as the initial guess; implemented by `pulse=[omega01_t; omega12_t]`.
- Lines 81-82: Run the optimisation, get normalised pulse; implemented by `fmaxnewton(spin_system,@grape_xy,pulse)`.

### Key state/data transformations

- Lines 16: computes `sys.magnet` using `sys.magnet=0`.
- Lines 19: computes `sys.isotopes` using `sys.isotopes={'T3'}`.
- Lines 22: computes `inter.modes.frqs` using `inter.modes.frqs={10e6}`.
- Lines 23: computes `inter.modes.anharms` using `inter.modes.anharms={-20e6}`.
- Lines 26: computes `bas.formalism` using `bas.formalism='zeeman-liouv'`.
- Lines 27: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 30: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 34: computes `H` using `H=hamiltonian(assume(spin_system,'cavity'))`.
- Lines 37: computes `pwr_levels` using `pwr_levels=2*pi*[30e6 35e6 40e6 45e6 50e6]`.
- Lines 40: computes `t` using `t=1e-9*linspace(-150,150,300)`.
- Lines 41: computes `ts` using `ts=-90e-9; sigma=45e-9`.
- Lines 42: computes `omega01` using `omega01=2*pi*43.4e6/max(pwr_levels)`.
- Lines 43: computes `omega12` using `omega12=2*pi*38.2e6/max(pwr_levels)`.
- Lines 44: computes `omega01_t` using `omega01_t=omega01*exp(-((t.^2)/(2*sigma^2)))`.
- Lines 45: computes `omega12_t` using `omega12_t=omega12*exp(-(((t-(ts/2)).^2)/(2*sigma^2)))`.
- Lines 48: computes `H01` using `H01=sparse([0 1 0; 1 0 0; 0 0 0]/2)`.
- Lines 49: computes `H12` using `H12=sparse([0 0 0; 0 0 1; 0 1 0]/2)`.
- Lines 52: computes `rho_init` using `rho_init=state(spin_system,'BL1',1)`.

## Implementation structure

- Basic single transmon system with Duffing model in-
- teractions, parameters and model from:
- Ensemble GRAPE optimisation with a distribution of
- control powers and a penalty on excess power.
- Calculation time: minutes
- Magnet field
- Particle specification
- Rotating frame ladder detunings
- Formalism and basis
- Spinach housekeeping
- Drift Hamiltonian from the declared interactions
- Pulse power ensemble

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `hamiltonian()`, `assume()`, `state()`, `hilb2liouv()`, `optimcon()`, `fmaxnewton()`.
