# examples/optimal_control/features_ampl.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/features_ampl.m`
- Signature: `features_ampl()`
- Total lines: 115

## Purpose

An illustration of amplitude profiling in a phase-modulated pulse optimisation. The amplitude profile is supplied by the user and the phase is optimised using LBFGS GRAPE algorithm with a penal- ty on the second derivative norm to encourage smoothness. In a set of 100 equspaced signals, the central 60 spins are set up for maximum excitation; there are no constraints on the dyna- mics of the 20 spins on either side of

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Set the magnetic field; implemented by `sys.magnet=14.1`.
- Lines 20-22: 100 non-interacting spins at equal intervals within the plus/minus 160 ppm range; implemented by `n_spins=100; sys.isotopes=cell(n_spins,1)`.
- Lines 26-28: Select a basis set -IK-2 keeps complete basis on each spin in this case, but ignores multi-spin orders; implemented by `bas.formalism='sphten-liouv'`.
- Lines 33-34: Run Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 37-38: Set up the initial state vector; implemented by `rho_init=state(spin_system,'Lz',21:80)`.
- Lines 41-42: Set up the target state vector; implemented by `rho_targ=state(spin_system,'Lx',21:80)`.
- Lines 45-46: Get the control operators; implemented by `Lx=operator(spin_system,'Lx','13C')`.
- Lines 49-50: Get the drift Hamiltonian; implemented by `D=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 52-53: Define control parameters; implemented by `control.isotopes={'13C'}`.
- Lines 67-69: Plotting options; implemented by `control.plotting={'xy_controls','phi_controls', 'robustness','spectrogram'}`.
- Lines 71-72: Initial guess for phases; implemented by `guess=pi/6*ones(1,250)`.
- Lines 74-75: Spinach housekeeping; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 77-78: Run LBFGS GRAPE pulse optimisation; implemented by `phi_profile=fmaxnewton(spin_system,@grape_phase,guess)`.
- Lines 80-81: Get Cartesian components of the pulse; implemented by `amp_profile=mean(control.pwr_levels)*control.amplitudes`.
- Lines 84-85: Simulate the optimised pulse; implemented by `rho_init=state(spin_system,'Lz','13C')`.
- Lines 89-90: Set acquisition parameters; implemented by `parameters.spins={'13C'}`.
- Lines 101-102: Simulate the free induction decay; implemented by `fid=liquid(spin_system,@acquire,parameters,'nmr')`.
- Lines 104-105: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'gauss',10}})`.

### Control flow inferred from the code

- Line 23: `for` loop over `n=1:n_spins, sys.isotopes{n}='13C'; end`.

### Key state/data transformations

- Lines 18: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 22: computes `n_spins` using `n_spins=100; sys.isotopes=cell(n_spins,1)`.
- Lines 24: computes `inter.zeeman.scalar` using `inter.zeeman.scalar=num2cell(linspace(-160,160,n_spins))`.
- Lines 28: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 29: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 30: computes `bas.space_level` using `bas.space_level=1`.
- Lines 31: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 34: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 38: computes `rho_init` using `rho_init=state(spin_system,'Lz',21:80)`.
- Lines 42: computes `rho_targ` using `rho_targ=state(spin_system,'Lx',21:80)`.
- Lines 46: computes `Lx` using `Lx=operator(spin_system,'Lx','13C')`.
- Lines 47: computes `Ly` using `Ly=operator(spin_system,'Ly','13C')`.
- Lines 50: computes `D` using `D=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 53: computes `control.isotopes` using `control.isotopes={'13C'}`.
- Lines 54: computes `control.channels` using `control.channels=[1; 1]`.
- Lines 55: computes `control.drifts` using `control.drifts={{D}}`.
- Lines 56: computes `control.operators` using `control.operators={Lx,Ly}`.
- Lines 57: computes `control.rho_init` using `control.rho_init={rho_init}`.

## Implementation structure

- An illustration of amplitude profiling in a phase-modulated pulse
- optimisation. The amplitude profile is supplied by the user and
- the phase is optimised using LBFGS GRAPE algorithm with a penal-
- ty on the second derivative norm to encourage smoothness.
- In a set of 100 equspaced signals, the central 60 spins are set
- up for maximum excitation; there are no constraints on the dyna-
- mics of the 20 spins on either side of the interval.
- Calculation time: minutes.
- Set the magnetic field
- 100 non-interacting spins at equal intervals
- within the plus/minus 160 ppm range
- Select a basis set -IK-2 keeps complete basis on each

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `num2cell()`, `create()`, `basis()`, `state()`, `operator()`, `hamiltonian()`, `assume()`, `optimcon()`, `fmaxnewton()`, `polar2cartesian()`, `shaped_pulse_xy()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
