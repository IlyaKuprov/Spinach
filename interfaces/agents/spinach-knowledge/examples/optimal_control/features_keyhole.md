# examples/optimal_control/features_keyhole.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/features_keyhole.m`
- Signature: `features_keyhole()`
- Total lines: 102

## Purpose

Optimal control pulse optimisation for state-to-state transfer across a scalar coupling in a hydrofluorocarbon fragment spin system. The start- ing state is Z-magnetisation on 1H, the destination state is Z-magneti- sation on 19F. There are six control channels and a distribution over control powers. A keyhole condition is specified: the system must be in a two-spin or- der at point 20 in the pulse time grid. LBFGS q

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 18-19: Magnetic field; implemented by `sys.magnet=9.4`.
- Lines 21-22: Spin system; implemented by `sys.isotopes={'1H','13C','19F'}`.
- Lines 24-25: Chemical shifts, ppm; implemented by `inter.zeeman.scalar={0.0, 0.0, 0.0}`.
- Lines 27-28: Scalar couplings, Hz (literature values); implemented by `inter.coupling.scalar=cell(3)`.
- Lines 32-33: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 36-37: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 40-41: Set up and normalise the initial state; implemented by `rho_init=state(spin_system,{'Lz'},{1})`.
- Lines 44-45: Set up and normalise the target state; implemented by `rho_targ=state(spin_system,{'Lz'},{3})`.
- Lines 48-49: Control operators; implemented by `LxH=operator(spin_system,'Lx','1H')`.
- Lines 56-57: Drift Hamiltonian; implemented by `D=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 59-60: Define control parameters; implemented by `control.isotopes={'1H','13C','19F'}`.
- Lines 73-74: Set the keyhole; implemented by `control.keyholes=cell(1,50)`.
- Lines 77-79: Plots during optimisation; implemented by `control.plotting={'correlation_order','local_each_spin', 'xy_controls','spectrogram'}`.
- Lines 81-82: Spinach housekeeping; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 84-85: Initial guess; implemented by `guess=randn(6,50)/3`.
- Lines 87-88: Run the optimal control procedure; implemented by `pulse=fmaxnewton(spin_system,@grape_xy,guess)`.
- Lines 90-91: Denormalise and format the waveform; implemented by `pulse=pulse*mean(control.pwr_levels)`.
- Lines 94-95: Run a test simulation using the optimal pulse; implemented by `report(spin_system,'running test simulation ')`.

### Key state/data transformations

- Lines 19: computes `sys.magnet` using `sys.magnet=9.4`.
- Lines 22: computes `sys.isotopes` using `sys.isotopes={'1H','13C','19F'}`.
- Lines 25: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.0, 0.0, 0.0}`.
- Lines 28: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(3)`.
- Lines 29: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=140`.
- Lines 30: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=-160`.
- Lines 33: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 34: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 37: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 41: computes `rho_init` using `rho_init=state(spin_system,{'Lz'},{1})`.
- Lines 45: computes `rho_targ` using `rho_targ=state(spin_system,{'Lz'},{3})`.
- Lines 49: computes `LxH` using `LxH=operator(spin_system,'Lx','1H')`.
- Lines 50: computes `LyH` using `LyH=operator(spin_system,'Ly','1H')`.
- Lines 51: computes `LxC` using `LxC=operator(spin_system,'Lx','13C')`.
- Lines 52: computes `LyC` using `LyC=operator(spin_system,'Ly','13C')`.
- Lines 53: computes `LxF` using `LxF=operator(spin_system,'Lx','19F')`.
- Lines 54: computes `LyF` using `LyF=operator(spin_system,'Ly','19F')`.
- Lines 57: computes `D` using `D=hamiltonian(assume(spin_system,'nmr'))`.

## Implementation structure

- Optimal control pulse optimisation for state-to-state transfer across a
- scalar coupling in a hydrofluorocarbon fragment spin system. The start-
- ing state is Z-magnetisation on 1H, the destination state is Z-magneti-
- sation on 19F. There are six control channels and a distribution over
- control powers.
- A keyhole condition is specified: the system must be in a two-spin or-
- der at point 20 in the pulse time grid. LBFGS quasi-Newton optimiser is
- used with a piecewise-linear control sequence.
- Calculation time: minutes.
- Magnetic field
- Spin system
- Chemical shifts, ppm

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `hamiltonian()`, `assume()`, `correlation()`, `optimcon()`, `fmaxnewton()`, `mat2cell()`, `report()`, `shaped_pulse_xy()`, `rho()`, `num2str()`.
