# examples/quantum_tech/transmon_frog.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/transmon_frog.m`
- Signature: `transmon_frog()`
- Total lines: 85

## Purpose

Basic implementation of a Frequency Robust Gate (FROG) for a single transmon, based on: Ensemble GRAPE optimisation with a distribution of control powers and excess amplitude penalty. Calculation time: minutes

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.
- The effective hardware model is a weakly anharmonic oscillator. Duffing nonlinearity breaks equal level spacing and allows qubit-like addressability within a truncated bosonic ladder.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 15-16: Magnet field; implemented by `sys.magnet=0`.
- Lines 18-19: Particle specification; implemented by `sys.isotopes={'T3'}`.
- Lines 21-22: Rotating frame transmon parameters; implemented by `inter.modes.frqs={0.5e6}`.
- Lines 25-26: Formalism and basis; implemented by `bas.formalism='zeeman-liouv'`.
- Lines 29-30: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 33-34: Drift Hamiltonian from the declared interactions; implemented by `H=hamiltonian(assume(spin_system,'cavity'))`.
- Lines 36-37: Build the intial control pulses (FROG); implemented by `nsteps=224; t_g=112e-9; t=linspace(0,t_g,nsteps)`.
- Lines 46-47: Quadrature control operators; implemented by `Cr=operator(spin_system,'C',1)`.
- Lines 51-52: Source and destination states; implemented by `psi_targ=[1; -1i; 0]/sqrt(2)`.
- Lines 58-59: Define control parameters; implemented by `control.isotopes={'T3'}`.
- Lines 72-73: Plots during optimisation; implemented by `control.plotting={'xy_controls','spectrogram','robustness'}`.
- Lines 75-76: Spinach housekeeping; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 78-79: Initial guess; implemented by `pulse=[omega_x; omega_y]`.
- Lines 81-82: Run the optimisation; implemented by `fmaxnewton(spin_system,@grape_xy,pulse)`.

### Control flow inferred from the code

- Line 41: `for` loop over `n=1:length(Fa)`.

### Key state/data transformations

- Lines 16: computes `sys.magnet` using `sys.magnet=0`.
- Lines 19: computes `sys.isotopes` using `sys.isotopes={'T3'}`.
- Lines 22: computes `inter.modes.frqs` using `inter.modes.frqs={0.5e6}`.
- Lines 23: computes `inter.modes.anharms` using `inter.modes.anharms={-295.1e6}`.
- Lines 26: computes `bas.formalism` using `bas.formalism='zeeman-liouv'`.
- Lines 27: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 30: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 34: computes `H` using `H=hamiltonian(assume(spin_system,'cavity'))`.
- Lines 37: computes `nsteps` using `nsteps=224; t_g=112e-9; t=linspace(0,t_g,nsteps)`.
- Lines 38: computes `Fa` using `Fa=[-0.6137 -0.0247 0.0742 0.0507 0.0149]`.
- Lines 39: computes `Fb` using `Fb=[-0.0106 0.0334 0.0579 0.0140 -0.0416]`.
- Lines 40: computes `omega_x` using `omega_x=zeros(1,nsteps); omega_y=zeros(1,nsteps)`.
- Lines 43: computes `omega_y` using `omega_y=omega_y+Fb(n)*sin((2*n*pi*t)/t_g)`.
- Lines 47: computes `Cr` using `Cr=operator(spin_system,'C',1)`.
- Lines 48: computes `An` using `An=operator(spin_system,'A',1)`.
- Lines 49: computes `Cx` using `Cx=(Cr+An)/2; Cy=1i*(Cr-An)/2`.
- Lines 52: computes `psi_targ` using `psi_targ=[1; -1i; 0]/sqrt(2)`.
- Lines 53: computes `rho_init` using `rho_init=state(spin_system,'BL1',1)`.

## Implementation structure

- Basic implementation of a Frequency Robust Gate
- (FROG) for a single transmon, based on:
- Ensemble GRAPE optimisation with a distribution
- of control powers and excess amplitude penalty.
- Calculation time: minutes
- Magnet field
- Particle specification
- Rotating frame transmon parameters
- Formalism and basis
- Spinach housekeeping
- Drift Hamiltonian from the declared interactions
- Build the intial control pulses (FROG)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `hamiltonian()`, `assume()`, `operator()`, `state()`, `hilb2liouv()`, `optimcon()`, `fmaxnewton()`.
