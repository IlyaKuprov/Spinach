# examples/optimal_control/features_diss_drift.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/features_diss_drift.m`
- Signature: `features_diss_drift()`
- Total lines: 135

## Purpose

Optimal control optimisation of a pulse performing magnetisa- tion transfer from H(N) to C(O) in a typical protein backbone spin system (literature data for shifts and couplings) with a range of pulse powers emulating B1 inhomogeneity and a range of offsets to account for imperfect transmitter placement. The dynamics includes dissipative terms in the drift generator: C(O) and N(H) are set to have rapid transverse rel

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

- Lines 24-25: Magnetic field; implemented by `sys.magnet=9.4`.
- Lines 27-28: Spin system; implemented by `sys.isotopes={'15N','1H','13C','13C','13C','15N'}`.
- Lines 31-32: Textbook chemical shifts, ppm; implemented by `inter.zeeman.scalar={119.79, 8.03, 57.32, 27.71, 177.25, 115.55}`.
- Lines 34-35: Scalar couplings, Hz (literature values); implemented by `inter.coupling.scalar=cell(6)`.
- Lines 43-44: Relaxation theory; implemented by `inter.relaxation={'t1_t2'}`.
- Lines 50-51: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 55-56: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 59-60: Set up and normalise the initial state; implemented by `rho_init=state(spin_system,{'Lz'},{2})`.
- Lines 63-64: Set up and normalise the target state; implemented by `rho_targ=state(spin_system,{'Lz'},{5})`.
- Lines 67-68: Control operators; implemented by `LxH=operator(spin_system,'Lx','1H')`.
- Lines 75-76: Offset operators; implemented by `LzH=operator(spin_system,'Lz','1H')`.
- Lines 80-81: Dissipative drift Liouvillian; implemented by `spin_system=assume(spin_system,'nmr')`.
- Lines 84-85: Put transmitters in the right place; implemented by `parameters.spins={'1H','13C','15N'}`.
- Lines 89-90: Define control parameters; implemented by `control.isotopes={'1H','13C','15N'}`.
- Lines 107-109: Control trajectory analysis plots; implemented by `control.plotting={'correlation_order','local_each_spin', 'amp_controls','spectrogram'}`.
- Lines 111-112: Spinach housekeeping; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 114-116: Start with pi/2 on 1H; end with pi/2 on 13C; implemented by `guess=randn(6,500)/10`.
- Lines 120-121: Run the optimisation; implemented by `pulse=fmaxnewton(spin_system,@grape_xy,guess)`.

### Key state/data transformations

- Lines 25: computes `sys.magnet` using `sys.magnet=9.4`.
- Lines 28: computes `sys.isotopes` using `sys.isotopes={'15N','1H','13C','13C','13C','15N'}`.
- Lines 29: computes `sys.labels` using `sys.labels={'N$_{(n)}$','H','CA','CB','C','N$_{(n+1)}$'}`.
- Lines 32: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={119.79, 8.03, 57.32, 27.71, 177.25, 115.55}`.
- Lines 35: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(6)`.
- Lines 36: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=-11`.
- Lines 37: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=140`.
- Lines 38: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4}=35`.
- Lines 39: computes `inter.coupling.scalar{3,5}` using `inter.coupling.scalar{3,5}=55`.
- Lines 40: computes `inter.coupling.scalar{3,6}` using `inter.coupling.scalar{3,6}=7`.
- Lines 41: computes `inter.coupling.scalar{5,6}` using `inter.coupling.scalar{5,6}=-15`.
- Lines 44: computes `inter.relaxation` using `inter.relaxation={'t1_t2'}`.
- Lines 45: computes `inter.r1_rates` using `inter.r1_rates={ 1.0 1.0 1.0 1.0 1.0 1.0}`.
- Lines 46: computes `inter.r2_rates` using `inter.r2_rates={50.0 1.0 1.0 1.0 100.0 50.0}`.
- Lines 47: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 48: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 51: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 52: computes `bas.approximation` using `bas.approximation='IK-0'`.

## Implementation structure

- Optimal control optimisation of a pulse performing magnetisa-
- tion transfer from H(N) to C(O) in a typical protein backbone
- spin system (literature data for shifts and couplings) with a
- range of pulse powers emulating B1 inhomogeneity and a range
- of offsets to account for imperfect transmitter placement.
- The dynamics includes dissipative terms in the drift generator:
- C(O) and N(H) are set to have rapid transverse relaxation. Four-
- spin correlation approximation is used, wherein five-spin and
- higher correlations are dropped from the basis set.
- The waveform is optimized with LBFGS-GRAPE algorithm with point-
- by-point variation and a penalty on the waveform amplitude.
- Calculation time: hours.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `assume()`, `hamiltonian()`, `relaxation()`, `frqoffset()`, `optimcon()`, `guess()`, `fmaxnewton()`, `mat2cell()`, `report()`, `shaped_pulse_xy()`, `rho()`, `num2str()`.
