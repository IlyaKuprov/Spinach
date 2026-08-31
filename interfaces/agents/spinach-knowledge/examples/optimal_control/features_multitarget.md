# examples/optimal_control/features_multitarget.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/features_multitarget.m`
- Signature: `features_multitarget()`
- Total lines: 110

## Purpose

An example of multi-target optimal control pulse design in the context of singlet state NMR spectroscopy. A pulse is designed that moves TT (carbon-triplet, proton-triplet) into SS (carbon-singlet, proton-singlet) and TS (carbon-triplet, proton-singlet) into ST (carbon-singlet, proton- triplet) simultaneously. The system is assumed to have a distribution in one of the J-couplings. Calculation time: hours.

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

- Lines 15-16: Magnetic field; implemented by `sys.magnet=14.1`.
- Lines 18-19: Isotopes; implemented by `sys.isotopes={'1H','13C','13C','1H'}`.
- Lines 21-22: Interactions; implemented by `inter.zeeman.scalar={0.0, 0.0, 0.0, 0.0}`.
- Lines 31-32: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 35-36: Run Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 39-40: Set up relevant operators; implemented by `unit_left=unit_oper(spin_system)`.
- Lines 52-53: Set up source and target states; implemented by `TT=CT_left*HT_left*unit_state(spin_system); TT=TT/norm(TT,2)`.
- Lines 59-60: Get the control operators; implemented by `LxH=operator(spin_system,'Lx','1H')`.
- Lines 65-66: Generate a distribution in J-couplings; implemented by `J=linspace(13,16,11); D=cell(numel(J),1)`.
- Lines 72-73: Define control parameters; implemented by `control.isotopes={'1H','13C'}`.
- Lines 86-88: Control trajectory analysis plots; implemented by `control.plotting={'correlation_order','local_each_spin', 'xy_controls','robustness','spectrogram'}`.
- Lines 90-91: Make control system structure; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 93-94: Run the optimisation from a random guess; implemented by `pulse=fmaxnewton(spin_system,@grape_xy,randn(4,275)/3)`.
- Lines 96-97: Denormalise and format the waveform; implemented by `pulse=pulse*control.pwr_levels`.
- Lines 100-101: Run a test simulation using the optimal pulse; implemented by `report(spin_system,'running test simulation ')`.

### Control flow inferred from the code

- Line 67: `for` loop over `m=1:numel(D)`.

### Key state/data transformations

- Lines 16: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 19: computes `sys.isotopes` using `sys.isotopes={'1H','13C','13C','1H'}`.
- Lines 22: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.0, 0.0, 0.0, 0.0}`.
- Lines 23: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(4)`.
- Lines 24: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=15.0`.
- Lines 25: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4}=15.0`.
- Lines 26: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=3.0`.
- Lines 27: computes `inter.coupling.scalar{2,4}` using `inter.coupling.scalar{2,4}=3.0`.
- Lines 28: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=150`.
- Lines 29: computes `inter.coupling.scalar{1,4}` using `inter.coupling.scalar{1,4}=8.0`.
- Lines 32: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 33: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 36: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 40: computes `unit_left` using `unit_left=unit_oper(spin_system)`.
- Lines 41: computes `HxHx_left` using `HxHx_left=operator(spin_system,{'Lx','Lx'},{1,4},'left')`.
- Lines 42: computes `HyHy_left` using `HyHy_left=operator(spin_system,{'Ly','Ly'},{1,4},'left')`.
- Lines 43: computes `HzHz_left` using `HzHz_left=operator(spin_system,{'Lz','Lz'},{1,4},'left')`.
- Lines 44: computes `CxCx_left` using `CxCx_left=operator(spin_system,{'Lx','Lx'},{2,3},'left')`.

## Implementation structure

- An example of multi-target optimal control pulse design in the context
- of singlet state NMR spectroscopy. A pulse is designed that moves TT
- (carbon-triplet, proton-triplet) into SS (carbon-singlet, proton-singlet)
- and TS (carbon-triplet, proton-singlet) into ST (carbon-singlet, proton-
- triplet) simultaneously. The system is assumed to have a distribution in
- one of the J-couplings.
- Calculation time: hours.
- Magnetic field
- Isotopes
- Interactions
- Basis set
- Run Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `unit_oper()`, `operator()`, `unit_state()`, `hamiltonian()`, `assume()`, `optimcon()`, `fmaxnewton()`, `mat2cell()`, `report()`, `shaped_pulse_xy()`, `cell2mat()`, `rho()`, `num2str()`.
