# examples/optimal_control/state_transfer_pro.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/state_transfer_pro.m`
- Signature: `state_transfer_pro()`
- Total lines: 122

## Purpose

Optimal control optimisation of a pulse performing magnetisa- tion transfer from H(N) to C(O) in a typical protein backbone spin system (literature data for shifts and couplings) with a range of pulse powers emulating B1 inhomogeneity and a range of offsets to account for imperfect transmitter placement. The waveform is optimized with LBFGS-GRAPE algorithm with po- int-by-point variation and a penalty on the pulse am

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

- Lines 19-20: Magnetic field; implemented by `sys.magnet=9.4`.
- Lines 22-23: Spin system; implemented by `sys.isotopes={'15N','1H','13C','13C','13C','15N'}`.
- Lines 26-27: Textbook chemical shifts, ppm; implemented by `inter.zeeman.scalar={119.79, 8.03, 57.32, 27.71, 177.25, 115.55}`.
- Lines 29-30: Scalar couplings, Hz (literature values); implemented by `inter.coupling.scalar=cell(6)`.
- Lines 38-39: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 43-44: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 47-48: Set up and normalise the initial state; implemented by `rho_init=state(spin_system,{'Lz'},{2})`.
- Lines 51-52: Set up and normalise the target state; implemented by `rho_targ=state(spin_system,{'Lz'},{5})`.
- Lines 55-56: Control operators; implemented by `LxH=operator(spin_system,'Lx','1H')`.
- Lines 63-64: Offset operators; implemented by `LzH=operator(spin_system,'Lz','1H')`.
- Lines 68-69: Drift Hamiltonian; implemented by `H=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 71-72: Put transmitters in the right place; implemented by `parameters.spins={'1H','13C','15N'}`.
- Lines 76-77: Define control parameters; implemented by `control.isotopes={'1H','13C','15N'}`.
- Lines 94-96: Control trajectory analysis plots; implemented by `control.plotting={'correlation_order','local_each_spin', 'amp_controls','spectrogram'}`.
- Lines 98-99: Spinach housekeeping; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 101-103: Start with pi/2 on 1H; end with pi/2 on 13C; implemented by `guess=randn(6,500)/10`.
- Lines 107-108: Run the optimisation; implemented by `pulse=fmaxnewton(spin_system,@grape_xy,guess)`.
- Lines 110-111: Apply power level scaling; implemented by `pulse=mean(control.pwr_levels)*pulse`.

### Key state/data transformations

- Lines 20: computes `sys.magnet` using `sys.magnet=9.4`.
- Lines 23: computes `sys.isotopes` using `sys.isotopes={'15N','1H','13C','13C','13C','15N'}`.
- Lines 24: computes `sys.labels` using `sys.labels={'N$_{(n)}$','H','CA','CB','C','N$_{(n+1)}$'}`.
- Lines 27: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={119.79, 8.03, 57.32, 27.71, 177.25, 115.55}`.
- Lines 30: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(6)`.
- Lines 31: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=-11`.
- Lines 32: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=140`.
- Lines 33: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4}=35`.
- Lines 34: computes `inter.coupling.scalar{3,5}` using `inter.coupling.scalar{3,5}=55`.
- Lines 35: computes `inter.coupling.scalar{3,6}` using `inter.coupling.scalar{3,6}=7`.
- Lines 36: computes `inter.coupling.scalar{5,6}` using `inter.coupling.scalar{5,6}=-15`.
- Lines 39: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 40: computes `bas.approximation` using `bas.approximation='IK-0'`.
- Lines 41: computes `bas.level` using `bas.level=4`.
- Lines 44: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 48: computes `rho_init` using `rho_init=state(spin_system,{'Lz'},{2})`.
- Lines 52: computes `rho_targ` using `rho_targ=state(spin_system,{'Lz'},{5})`.
- Lines 56: computes `LxH` using `LxH=operator(spin_system,'Lx','1H')`.

## Implementation structure

- Optimal control optimisation of a pulse performing magnetisa-
- tion transfer from H(N) to C(O) in a typical protein backbone
- spin system (literature data for shifts and couplings) with a
- range of pulse powers emulating B1 inhomogeneity and a range
- of offsets to account for imperfect transmitter placement.
- The waveform is optimized with LBFGS-GRAPE algorithm with po-
- int-by-point variation and a penalty on the pulse amplitude.
- Calculation time: hours.
- Magnetic field
- Spin system
- Textbook chemical shifts, ppm
- Scalar couplings, Hz (literature values)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `hamiltonian()`, `assume()`, `frqoffset()`, `optimcon()`, `guess()`, `fmaxnewton()`, `mat2cell()`, `report()`, `shaped_pulse_xy()`, `rho()`, `num2str()`.
