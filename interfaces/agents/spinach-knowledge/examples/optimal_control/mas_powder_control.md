# examples/optimal_control/mas_powder_control.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/mas_powder_control.m`
- Signature: `mas_powder_control()`
- Total lines: 106

## Purpose

Optimal control pulse starting with Lz and populating the Ly state on 87Rb in a quadrupolar rubidium system under magic angle spinning. A phase-modulated pulse is produced. Calculation time: hours.

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: System specification; implemented by `sys.magnet=9.413`.
- Lines 16-17: Quadrupolar coupling; implemented by `inter.coupling.matrix{1,1}=eeqq2nqi(1.68e6,0.2,3/2,[0 0 0])`.
- Lines 19-20: Basis set and formalism; implemented by `bas.formalism='sphten-liouv'`.
- Lines 23-24: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 27-28: MAS experiment parameters; implemented by `parameters.axis=[sqrt(2/3) 0 sqrt(1/3)]`.
- Lines 36-37: Drift Liouvillians and classical subspace dimension for the ensemble; implemented by `[control.drifts,spc_dim]=drifts(spin_system,@singlerot,parameters,'qnmr')`.
- Lines 39-40: Initial state -Lz; implemented by `rho_init=state(spin_system,'Lz','87Rb')`.
- Lines 46-47: Target state -Ly; implemented by `rho_targ=state(spin_system,'Ly','87Rb')`.
- Lines 52-53: Control operators; implemented by `Lx=operator(spin_system,'Lx','87Rb')`.
- Lines 59-60: Offset operator; implemented by `Lz=operator(spin_system,'Lz','87Rb')`.
- Lines 63-64: Control parameters; implemented by `control.isotopes={'87Rb'}`.
- Lines 75-77: Plotting options; implemented by `control.plotting={'phi_controls','robustness', 'xy_controls','spectrogram'}`.
- Lines 79-80: Initial guess; implemented by `guess=pi/2*ones(1,100)`.
- Lines 82-83: Spinach housekeeping; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 85-86: Run the optimisation; implemented by `pulse=fmaxnewton(spin_system,@grape_phase,guess)`.
- Lines 88-89: Construct Cartesian waveform; implemented by `pulse_x=mean(control.pwr_levels)*cos(pulse)`.
- Lines 92-93: Average over the systems; implemented by `parfor n=1:numel(control.drifts)`.
- Lines 95-97: Run the pulse; implemented by `rho=shaped_pulse_xy(spin_system,control.drifts{n}{1},control.operators, {pulse_x,pulse_y},control.pulse_dt,rho_init,'expv-pwc')`.

### Control flow inferred from the code

- Line 93: `parfor` loop over `n=1:numel(control.drifts)`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=9.413`.
- Lines 14: computes `sys.isotopes` using `sys.isotopes={'87Rb'}`.
- Lines 17: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=eeqq2nqi(1.68e6,0.2,3/2,[0 0 0])`.
- Lines 20: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 21: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 24: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 28: computes `parameters.axis` using `parameters.axis=[sqrt(2/3) 0 sqrt(1/3)]`.
- Lines 29: computes `parameters.rate` using `parameters.rate=-20e3`.
- Lines 30: computes `parameters.grid` using `parameters.grid='rep_2ang_100pts_sph'`.
- Lines 31: computes `parameters.max_rank` using `parameters.max_rank=8`.
- Lines 32: computes `parameters.spins` using `parameters.spins={'87Rb'}`.
- Lines 33: computes `parameters.rframes` using `parameters.rframes={{'87Rb',3}}`.
- Lines 34: computes `parameters.verbose` using `parameters.verbose=0`.
- Lines 37: computes `[control.drifts,spc_dim]` using `[control.drifts,spc_dim]=drifts(spin_system,@singlerot,parameters,'qnmr')`.
- Lines 40: computes `rho_init` using `rho_init=state(spin_system,'Lz','87Rb')`.
- Lines 41: computes `space_part` using `space_part=ones(spc_dim,1)`.
- Lines 44: computes `control.rho_init` using `control.rho_init={rho_init}`.
- Lines 47: computes `rho_targ` using `rho_targ=state(spin_system,'Ly','87Rb')`.

## Implementation structure

- Optimal control pulse starting with Lz and populating the
- Ly state on 87Rb in a quadrupolar rubidium system under
- magic angle spinning. A phase-modulated pulse is produced.
- Calculation time: hours.
- System specification
- Quadrupolar coupling
- Basis set and formalism
- Spinach housekeeping
- MAS experiment parameters
- Drift Liouvillians and classical subspace dimension for the ensemble
- Initial state -Lz
- Target state -Ly

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `drifts()`, `state()`, `operator()`, `speye()`, `optimcon()`, `fmaxnewton()`, `shaped_pulse_xy()`, `fid()`, `num2str()`.
