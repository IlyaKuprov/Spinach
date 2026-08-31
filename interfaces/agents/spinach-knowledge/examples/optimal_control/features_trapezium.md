# examples/optimal_control/features_trapezium.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/features_trapezium.m`
- Signature: `features_trapezium()`
- Total lines: 113

## Purpose

Optimal control pulse optimisation for state-to-state transfer across a scalar coupling in a hydrofluorocarbon fragment spin system. The start- ing state is Z-magnetisation on 1H, the destination state is Z-magneti- sation on 19F. There are six control channels. The waveform is treated as piecewise-linear using derivatives of a Lie group product quadrature, as described in: Calculation time: minutes.

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 19-20: Magnetic field; implemented by `sys.magnet=9.4`.
- Lines 22-23: Spin system; implemented by `sys.isotopes={'1H','13C','19F'}`.
- Lines 25-26: Chemical shifts, ppm; implemented by `inter.zeeman.scalar={0.0, 0.0, 0.0}`.
- Lines 28-29: Scalar couplings, Hz (literature values); implemented by `inter.coupling.scalar=cell(3)`.
- Lines 33-34: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 37-38: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 41-42: Set up and normalise the initial state; implemented by `rho_init=state(spin_system,{'Lz'},{1})`.
- Lines 45-46: Set up and normalise the target state; implemented by `rho_targ=state(spin_system,{'Lz'},{3})`.
- Lines 49-50: Control operators; implemented by `LxH=operator(spin_system,'Lx','1H')`.
- Lines 57-58: Drift Hamiltonian; implemented by `D=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 60-61: Define control parameters; implemented by `control.isotopes={'1H','13C','19F'}`.
- Lines 75-77: Plots during optimisation; implemented by `control.plotting={'correlation_order','local_each_spin', 'xy_controls','spectrogram'}`.
- Lines 79-80: Spinach housekeeping; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 82-83: Initial guess; implemented by `guess=randn(6,51)/3`.
- Lines 85-86: Run the optimal control procedure; implemented by `pulse=fmaxnewton(spin_system,@grape_xy,guess)`.
- Lines 88-89: Denormalise the waveform; implemented by `pulse=pulse*mean(control.pwr_levels)`.
- Lines 91-92: Run a test simulation using the optimal pulse; implemented by `report(spin_system,'running test simulation ')`.
- Lines 96-97: Build generators; implemented by `G_L=D; G_R=D`.

### Control flow inferred from the code

- Line 94: `for` loop over `n=1:numel(control.pulse_dt)`.
- Line 98: `for` loop over `k=1:size(pulse,1)`.

### Key state/data transformations

- Lines 20: computes `sys.magnet` using `sys.magnet=9.4`.
- Lines 23: computes `sys.isotopes` using `sys.isotopes={'1H','13C','19F'}`.
- Lines 26: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.0, 0.0, 0.0}`.
- Lines 29: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(3)`.
- Lines 30: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=140`.
- Lines 31: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=-160`.
- Lines 34: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 35: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 38: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 42: computes `rho_init` using `rho_init=state(spin_system,{'Lz'},{1})`.
- Lines 46: computes `rho_targ` using `rho_targ=state(spin_system,{'Lz'},{3})`.
- Lines 50: computes `LxH` using `LxH=operator(spin_system,'Lx','1H')`.
- Lines 51: computes `LyH` using `LyH=operator(spin_system,'Ly','1H')`.
- Lines 52: computes `LxC` using `LxC=operator(spin_system,'Lx','13C')`.
- Lines 53: computes `LyC` using `LyC=operator(spin_system,'Ly','13C')`.
- Lines 54: computes `LxF` using `LxF=operator(spin_system,'Lx','19F')`.
- Lines 55: computes `LyF` using `LyF=operator(spin_system,'Ly','19F')`.
- Lines 58: computes `D` using `D=hamiltonian(assume(spin_system,'nmr'))`.

## Implementation structure

- Optimal control pulse optimisation for state-to-state transfer across a
- scalar coupling in a hydrofluorocarbon fragment spin system. The start-
- ing state is Z-magnetisation on 1H, the destination state is Z-magneti-
- sation on 19F. There are six control channels.
- The waveform is treated as piecewise-linear using derivatives of a Lie
- group product quadrature, as described in:
- Calculation time: minutes.
- Magnetic field
- Spin system
- Chemical shifts, ppm
- Scalar couplings, Hz (literature values)
- Basis set

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `hamiltonian()`, `assume()`, `optimcon()`, `fmaxnewton()`, `report()`, `pulse()`, `step()`, `rho()`, `num2str()`.
