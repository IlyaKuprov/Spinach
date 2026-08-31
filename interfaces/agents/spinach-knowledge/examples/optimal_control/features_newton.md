# examples/optimal_control/features_newton.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/features_newton.m`
- Signature: `features_newton()`
- Total lines: 106

## Purpose

Optimal control pulse optimisation for state-to-state transfer across two scalar couplings in a hydrofluorocarbon fragment spin system. The starting state is Lz on 1H, the destination state is Lz on 19F. There are six control channels, the pulse is designed to be stable with res- pect to proton transmitter offset and pulse nutation frequency drop. The optimisation uses Newton-Raphson GRAPE algorithm described in: wit

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 21-22: Magnetic field; implemented by `sys.magnet=9.4`.
- Lines 24-25: Spin system; implemented by `sys.isotopes={'1H','13C','19F'}`.
- Lines 27-28: Chemical shifts, ppm; implemented by `inter.zeeman.scalar={0.0, 0.0, 0.0}`.
- Lines 30-31: Scalar couplings, Hz (literature values); implemented by `inter.coupling.scalar=cell(3)`.
- Lines 35-36: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 39-40: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 43-44: Set up and normalise the initial state; implemented by `rho_init=state(spin_system,{'Lz'},{1})`.
- Lines 47-48: Set up and normalise the target state; implemented by `rho_targ=state(spin_system,{'Lz'},{3})`.
- Lines 51-52: Control operators; implemented by `LxH=operator(spin_system,'Lx','1H')`.
- Lines 59-60: Offset operator; implemented by `LzH=operator(spin_system,'Lz','1H')`.
- Lines 62-63: Drift Hamiltonian; implemented by `D=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 65-66: Define control parameters; implemented by `control.isotopes={'1H','13C','19F'}`.
- Lines 81-83: Plots during optimisation; implemented by `control.plotting={'correlation_order','local_each_spin', 'xy_controls','spectrogram'}`.
- Lines 85-86: Spinach housekeeping; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 88-89: Initial guess; implemented by `guess=randn(6,100)/10`.
- Lines 91-92: Run the optimal control procedure; implemented by `pulse=fmaxnewton(spin_system,@grape_xy,guess)`.
- Lines 94-95: Denormalise and format the waveform; implemented by `pulse=pulse*mean(control.pwr_levels)`.
- Lines 98-99: Run a test simulation using the optimal pulse; implemented by `report(spin_system,'running test simulation ')`.

### Key state/data transformations

- Lines 22: computes `sys.magnet` using `sys.magnet=9.4`.
- Lines 25: computes `sys.isotopes` using `sys.isotopes={'1H','13C','19F'}`.
- Lines 28: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.0, 0.0, 0.0}`.
- Lines 31: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(3)`.
- Lines 32: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=140`.
- Lines 33: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=-160`.
- Lines 36: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 37: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 40: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 44: computes `rho_init` using `rho_init=state(spin_system,{'Lz'},{1})`.
- Lines 48: computes `rho_targ` using `rho_targ=state(spin_system,{'Lz'},{3})`.
- Lines 52: computes `LxH` using `LxH=operator(spin_system,'Lx','1H')`.
- Lines 53: computes `LyH` using `LyH=operator(spin_system,'Ly','1H')`.
- Lines 54: computes `LxC` using `LxC=operator(spin_system,'Lx','13C')`.
- Lines 55: computes `LyC` using `LyC=operator(spin_system,'Ly','13C')`.
- Lines 56: computes `LxF` using `LxF=operator(spin_system,'Lx','19F')`.
- Lines 57: computes `LyF` using `LyF=operator(spin_system,'Ly','19F')`.
- Lines 60: computes `LzH` using `LzH=operator(spin_system,'Lz','1H')`.

## Implementation structure

- Optimal control pulse optimisation for state-to-state transfer across
- two scalar couplings in a hydrofluorocarbon fragment spin system. The
- starting state is Lz on 1H, the destination state is Lz on 19F. There
- are six control channels, the pulse is designed to be stable with res-
- pect to proton transmitter offset and pulse nutation frequency drop.
- The optimisation uses Newton-Raphson GRAPE algorithm described in:
- with point-by-point variation and a penalty on the waveform exceeding
- a user-specified power threshold. The initial guess is a random pulse.
- Calculation time: minutes.
- Magnetic field
- Spin system
- Chemical shifts, ppm

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `hamiltonian()`, `assume()`, `optimcon()`, `fmaxnewton()`, `mat2cell()`, `report()`, `shaped_pulse_xy()`, `rho()`, `num2str()`.
