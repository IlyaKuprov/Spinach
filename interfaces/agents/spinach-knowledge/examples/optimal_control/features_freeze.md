# examples/optimal_control/features_freeze.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/features_freeze.m`
- Signature: `features_freeze()`
- Total lines: 114

## Purpose

Optimal control pulse optimisation for state-to-state transfer across a scalar coupling in a hydrofluorocarbon fragment spin system. The start- ing state is Z-magnetisation on 1H, the destination state is Z-magneti- sation on 19F. There are six control channels. A freeze condition is specified -there are two periods in the control sequence that the optimisation is not allowed to touch. The waveform is optimised with 

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 24-25: Magnetic field; implemented by `sys.magnet=9.4`.
- Lines 27-28: Spin system; implemented by `sys.isotopes={'1H','13C','19F'}`.
- Lines 30-31: Chemical shifts, ppm; implemented by `inter.zeeman.scalar={0.0, 0.0, 0.0}`.
- Lines 33-34: Scalar couplings, Hz (literature values); implemented by `inter.coupling.scalar=cell(3)`.
- Lines 38-39: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 42-43: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 46-47: Set up and normalise the initial state; implemented by `rho_init=state(spin_system,{'Lz'},{1})`.
- Lines 50-51: Set up and normalise the target state; implemented by `rho_targ=state(spin_system,{'Lz'},{3})`.
- Lines 54-55: Control operators; implemented by `LxH=operator(spin_system,'Lx','1H')`.
- Lines 62-63: Offset operator; implemented by `LzH=operator(spin_system,'Lz','1H')`.
- Lines 65-66: Drift Hamiltonian; implemented by `D=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 68-69: Define control parameters; implemented by `control.isotopes={'1H','13C','19F'}`.
- Lines 84-85: Initial guess; implemented by `guess=randn(6,100)/10`.
- Lines 88-89: Freeze region; implemented by `control.freeze=false(6,100)`.
- Lines 92-94: Plots during optimisation; implemented by `control.plotting={'correlation_order','local_each_spin', 'xy_controls','spectrogram'}`.
- Lines 96-97: Spinach housekeeping; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 99-100: Run the optimal control procedure; implemented by `pulse=fmaxnewton(spin_system,@grape_xy,guess)`.
- Lines 102-103: Denormalise and format the waveform; implemented by `pulse=pulse*mean(control.pwr_levels)`.

### Key state/data transformations

- Lines 25: computes `sys.magnet` using `sys.magnet=9.4`.
- Lines 28: computes `sys.isotopes` using `sys.isotopes={'1H','13C','19F'}`.
- Lines 31: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.0, 0.0, 0.0}`.
- Lines 34: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(3)`.
- Lines 35: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=140`.
- Lines 36: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=-160`.
- Lines 39: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 40: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 43: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 47: computes `rho_init` using `rho_init=state(spin_system,{'Lz'},{1})`.
- Lines 51: computes `rho_targ` using `rho_targ=state(spin_system,{'Lz'},{3})`.
- Lines 55: computes `LxH` using `LxH=operator(spin_system,'Lx','1H')`.
- Lines 56: computes `LyH` using `LyH=operator(spin_system,'Ly','1H')`.
- Lines 57: computes `LxC` using `LxC=operator(spin_system,'Lx','13C')`.
- Lines 58: computes `LyC` using `LyC=operator(spin_system,'Ly','13C')`.
- Lines 59: computes `LxF` using `LxF=operator(spin_system,'Lx','19F')`.
- Lines 60: computes `LyF` using `LyF=operator(spin_system,'Ly','19F')`.
- Lines 63: computes `LzH` using `LzH=operator(spin_system,'Lz','1H')`.

## Implementation structure

- Optimal control pulse optimisation for state-to-state transfer across a
- scalar coupling in a hydrofluorocarbon fragment spin system. The start-
- ing state is Z-magnetisation on 1H, the destination state is Z-magneti-
- sation on 19F. There are six control channels.
- A freeze condition is specified -there are two periods in the control
- sequence that the optimisation is not allowed to touch.
- The waveform is optimised with the Newton-Raphson GRAPE algorithm desc-
- ribed in
- with point-by-point variation and a penalty on the waveform exceeding
- a user-specified power threshold. The initial guess is a random pulse;
- the optimisation typically achieves a fidelity of 0.999999.
- Calculation time: minutes.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `hamiltonian()`, `assume()`, `guess()`, `false()`, `optimcon()`, `fmaxnewton()`, `mat2cell()`, `report()`, `shaped_pulse_xy()`, `rho()`, `num2str()`.
