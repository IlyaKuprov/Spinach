# examples/quantum_tech/transmon_transfer.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/transmon_transfer.m`
- Signature: `transmon_transfer()`
- Total lines: 95

## Purpose

Basic two-transmon system with Duffing model interacti- ons and a flip-flop coupling; coherence transfer from transmon 1 to transmon 2. GRAPE optimisation with a distribution control powers and transmon offsets with a penalty on excess power. Calculation time: minutes.

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.
- The effective hardware model is a weakly anharmonic oscillator. Duffing nonlinearity breaks equal level spacing and allows qubit-like addressability within a truncated bosonic ladder.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: Magnet field; implemented by `sys.magnet=0`.
- Lines 17-18: Particle specification; implemented by `sys.isotopes={'T3','T5'}`.
- Lines 20-21: Rotating frame transmon parameters; implemented by `inter.modes.frqs={100e6 -200e6}`.
- Lines 26-27: Formalism and basis; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 30-31: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 34-35: Drift Hamiltonian from the declared interactions; implemented by `H=hamiltonian(assume(spin_system,'cavity'))`.
- Lines 37-38: Get elementary operators; implemented by `CrA=operator(spin_system,'C',1)`.
- Lines 43-44: Build control operators; implemented by `C_A=(CrA+AnA)/2; C_B=(CrB+AnB)/2`.
- Lines 46-47: Build offset operators; implemented by `O_A=operator(spin_system,'N',1)`.
- Lines 50-52: Build source and destination states; implemented by `rho_init=state(spin_system,{'C','BL1'},{1 2})+ state(spin_system,{'A','BL1'},{1 2})`.
- Lines 58-59: Clean up rounding errors; implemented by `rho_init=(rho_init+rho_init')/2`.
- Lines 62-63: Unit fidelity is Sorensen bound; implemented by `rho_targ=rho_targ/sorensen(rho_init,rho_targ)`.
- Lines 65-66: Define control parameters; implemented by `control.isotopes={'T3','T5'}`.
- Lines 82-83: Plots during optimisation; implemented by `control.plotting={'xy_controls','spectrogram','robustness'}`.
- Lines 85-86: Spinach housekeeping; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 88-89: Initial guess -random; implemented by `pulse=(1/10)*randn(2,200)`.
- Lines 91-92: Run the optimisation, get normalised pulse; implemented by `fmaxnewton(spin_system,@grape_xy,pulse)`.

### Key state/data transformations

- Lines 15: computes `sys.magnet` using `sys.magnet=0`.
- Lines 18: computes `sys.isotopes` using `sys.isotopes={'T3','T5'}`.
- Lines 21: computes `inter.modes.frqs` using `inter.modes.frqs={100e6 -200e6}`.
- Lines 22: computes `inter.modes.anharms` using `inter.modes.anharms={-10e6 -20e6}`.
- Lines 23: computes `inter.modes.exchange` using `inter.modes.exchange=cell(2,2)`.
- Lines 24: computes `inter.modes.exchange{1,2}` using `inter.modes.exchange{1,2}=50e6`.
- Lines 27: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 28: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 31: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 35: computes `H` using `H=hamiltonian(assume(spin_system,'cavity'))`.
- Lines 38: computes `CrA` using `CrA=operator(spin_system,'C',1)`.
- Lines 39: computes `AnA` using `AnA=operator(spin_system,'A',1)`.
- Lines 40: computes `CrB` using `CrB=operator(spin_system,'C',2)`.
- Lines 41: computes `AnB` using `AnB=operator(spin_system,'A',2)`.
- Lines 44: computes `C_A` using `C_A=(CrA+AnA)/2; C_B=(CrB+AnB)/2`.
- Lines 47: computes `O_A` using `O_A=operator(spin_system,'N',1)`.
- Lines 48: computes `O_B` using `O_B=operator(spin_system,'N',2)`.
- Lines 51-52: computes `rho_init` using `rho_init=state(spin_system,{'C','BL1'},{1 2})+ state(spin_system,{'A','BL1'},{1 2})`.

## Implementation structure

- Basic two-transmon system with Duffing model interacti-
- ons and a flip-flop coupling; coherence transfer from
- transmon 1 to transmon 2. GRAPE optimisation with a
- distribution control powers and transmon offsets with
- a penalty on excess power.
- Calculation time: minutes.
- Magnet field
- Particle specification
- Rotating frame transmon parameters
- Formalism and basis
- Spinach housekeeping
- Drift Hamiltonian from the declared interactions

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `hamiltonian()`, `assume()`, `operator()`, `state()`, `sorensen()`, `optimcon()`, `fmaxnewton()`.
