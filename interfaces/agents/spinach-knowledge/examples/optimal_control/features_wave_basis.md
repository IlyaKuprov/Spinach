# examples/optimal_control/features_wave_basis.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/features_wave_basis.m`
- Signature: `features_wave_basis()`
- Total lines: 116

## Purpose

An illustration of basis set coefficient optimisation for a pul- se optimised as a linear combination of user-specified vectors. The optimisation is done using LBFGS GRAPE algorithm. In a set of 100 equspaced signals, the central 60 spins are set up for maximum excitation; there are no constraints on the dyna- mics of the 20 spins on either side of the interval. Calculation time: minutes.

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

- Lines 18-19: Set the magnetic field; implemented by `sys.magnet=14.1`.
- Lines 21-23: 100 non-interacting spins at equal intervals within the plus/minus 160 ppm range; implemented by `n_spins=100; sys.isotopes=cell(n_spins,1)`.
- Lines 27-29: Select a basis set -IK-2 keeps complete basis on each spin in this case, but ignores multi-spin orders; implemented by `bas.formalism='sphten-liouv'`.
- Lines 34-35: Run Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 38-39: Set up the initial state vector; implemented by `rho_init=state(spin_system,'Lz',21:80)`.
- Lines 42-43: Set up the target state vector; implemented by `rho_targ=state(spin_system,'Lx',21:80)`.
- Lines 46-47: Get the control operators; implemented by `Cx=operator(spin_system,'Lx','13C')`.
- Lines 50-51: Get the drift Hamiltonian; implemented by `D=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 53-54: Define control parameters; implemented by `control.isotopes={'13C'}`.
- Lines 68-70: Plotting options; implemented by `control.plotting={'xy_controls','amp_controls', 'robustness','spectrogram'}`.
- Lines 72-73: Initial guess; implemented by `guess=randn(2,20)/20`.
- Lines 75-76: Spinach housekeeping; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 78-79: Run LBFGS GRAPE pulse optimisation; implemented by `basis_coeffs=fmaxnewton(spin_system,@grape_xy,guess)`.
- Lines 81-82: Reassemble time-domain control sequence; implemented by `pulse=mean(control.pwr_levels)*basis_coeffs*control.basis`.
- Lines 85-86: Simulate the optimised pulse; implemented by `rho_init=state(spin_system,'Lz','13C')`.
- Lines 90-91: Set acquisition parameters; implemented by `parameters.spins={'13C'}`.
- Lines 102-103: Simulate the free induction decay; implemented by `fid=liquid(spin_system,@acquire,parameters,'nmr')`.
- Lines 105-106: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'gauss',10}})`.

### Control flow inferred from the code

- Line 24: `for` loop over `n=1:n_spins, sys.isotopes{n}='13C'; end`.

### Key state/data transformations

- Lines 19: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 23: computes `n_spins` using `n_spins=100; sys.isotopes=cell(n_spins,1)`.
- Lines 25: computes `inter.zeeman.scalar` using `inter.zeeman.scalar=num2cell(linspace(-160,160,n_spins))`.
- Lines 29: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 30: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 31: computes `bas.space_level` using `bas.space_level=1`.
- Lines 32: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 35: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 39: computes `rho_init` using `rho_init=state(spin_system,'Lz',21:80)`.
- Lines 43: computes `rho_targ` using `rho_targ=state(spin_system,'Lx',21:80)`.
- Lines 47: computes `Cx` using `Cx=operator(spin_system,'Lx','13C')`.
- Lines 48: computes `Cy` using `Cy=operator(spin_system,'Ly','13C')`.
- Lines 51: computes `D` using `D=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 54: computes `control.isotopes` using `control.isotopes={'13C'}`.
- Lines 55: computes `control.channels` using `control.channels=[1; 1]`.
- Lines 56: computes `control.drifts` using `control.drifts={{D}}`.
- Lines 57: computes `control.operators` using `control.operators={Cx,Cy}`.
- Lines 58: computes `control.rho_init` using `control.rho_init={rho_init}`.

## Implementation structure

- An illustration of basis set coefficient optimisation for a pul-
- se optimised as a linear combination of user-specified vectors.
- The optimisation is done using LBFGS GRAPE algorithm.
- In a set of 100 equspaced signals, the central 60 spins are set
- up for maximum excitation; there are no constraints on the dyna-
- mics of the 20 spins on either side of the interval.
- Calculation time: minutes.
- Set the magnetic field
- 100 non-interacting spins at equal intervals
- within the plus/minus 160 ppm range
- Select a basis set -IK-2 keeps complete basis on each
- spin in this case, but ignores multi-spin orders

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `num2cell()`, `create()`, `basis()`, `state()`, `operator()`, `hamiltonian()`, `assume()`, `wave_basis()`, `optimcon()`, `fmaxnewton()`, `mat2cell()`, `shaped_pulse_xy()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`.
