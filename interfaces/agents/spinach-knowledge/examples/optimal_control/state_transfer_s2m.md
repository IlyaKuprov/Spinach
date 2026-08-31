# examples/optimal_control/state_transfer_s2m.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/state_transfer_s2m.m`
- Signature: `state_transfer_s2m()`
- Total lines: 100

## Purpose

A transfer of coherence from a two-proton singlet state to a nearby carbon in a setting typically encountered in parahydrogenation expe- riments. LBFGS-GRAPE algorithm is used as described in Terminal fidelity is 50% in this case. Calculation time: minutes.

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

- Lines 16-17: Magnetic field; implemented by `sys.magnet=14.1`.
- Lines 19-20: Isotopes; implemented by `sys.isotopes={'1H','1H','13C','13C'}`.
- Lines 22-23: Interactions; implemented by `inter.zeeman.scalar={1.5, 2.0, 30.0, 40.0}`.
- Lines 30-31: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 34-35: Run Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 38-39: Set up and normalise the initial state; implemented by `rho_init=singlet(spin_system,1,2)`.
- Lines 42-43: Set up and normalise the target state; implemented by `rho_targ=state(spin_system,{'Lz'},{4})`.
- Lines 46-47: Get the control operators; implemented by `LxH=operator(spin_system,'Lx','1H')`.
- Lines 52-53: Drift Hamiltonian; implemented by `H=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 55-56: Transmitter offsets; implemented by `parameters.spins={'1H','13C'}`.
- Lines 60-61: Define control parameters; implemented by `control.isotopes={'1H','13C'}`.
- Lines 74-77: Plots during optimisation; implemented by `control.plotting={'correlation_order','coherence_order', 'xy_controls','local_each_spin', 'spectrogram'}`.
- Lines 79-80: Spinach housekeeping; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 82-83: Initial guess; implemented by `guess=randn(4,100)/10`.
- Lines 85-86: Run the optimisation, get normalised pulse; implemented by `pulse=fmaxnewton(spin_system,@grape_xy,guess)`.
- Lines 88-89: Apply power level scaling; implemented by `pulse=mean(control.pwr_levels)*pulse`.
- Lines 92-93: Run a test simulation using the optimal pulse; implemented by `report(spin_system,'running test simulation ')`.

### Key state/data transformations

- Lines 17: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 20: computes `sys.isotopes` using `sys.isotopes={'1H','1H','13C','13C'}`.
- Lines 23: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={1.5, 2.0, 30.0, 40.0}`.
- Lines 24: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(4)`.
- Lines 25: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=7.0`.
- Lines 26: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=150`.
- Lines 27: computes `inter.coupling.scalar{2,4}` using `inter.coupling.scalar{2,4}=150`.
- Lines 28: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4}=50`.
- Lines 31: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 32: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 35: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 39: computes `rho_init` using `rho_init=singlet(spin_system,1,2)`.
- Lines 43: computes `rho_targ` using `rho_targ=state(spin_system,{'Lz'},{4})`.
- Lines 47: computes `LxH` using `LxH=operator(spin_system,'Lx','1H')`.
- Lines 48: computes `LyH` using `LyH=operator(spin_system,'Ly','1H')`.
- Lines 49: computes `LxC` using `LxC=operator(spin_system,'Lx','13C')`.
- Lines 50: computes `LyC` using `LyC=operator(spin_system,'Ly','13C')`.
- Lines 53: computes `H` using `H=hamiltonian(assume(spin_system,'nmr'))`.

## Implementation structure

- A transfer of coherence from a two-proton singlet state to a nearby
- carbon in a setting typically encountered in parahydrogenation expe-
- riments. LBFGS-GRAPE algorithm is used as described in
- Terminal fidelity is 50% in this case.
- Calculation time: minutes.
- Magnetic field
- Isotopes
- Interactions
- Basis set
- Run Spinach housekeeping
- Set up and normalise the initial state
- Set up and normalise the target state

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `singlet()`, `state()`, `operator()`, `hamiltonian()`, `assume()`, `frqoffset()`, `optimcon()`, `fmaxnewton()`, `mat2cell()`, `report()`, `shaped_pulse_xy()`, `rho()`, `num2str()`.
