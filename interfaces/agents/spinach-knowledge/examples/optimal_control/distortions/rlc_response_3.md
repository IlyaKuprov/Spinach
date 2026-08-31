# examples/optimal_control/distortions/rlc_response_3.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/distortions/rlc_response_3.m`
- Signature: `rlc_response_3()`
- Total lines: 102

## Purpose

Probe circuit response effect on the accuracy of the deu- terium pre-phasing pulse designed to set deuterium magne- tisation in a -CD3 group of alanine up for rephasing 100 microseconds after the pulse is finished. The system is assumed to be a powder (100 orientations) with a B1 distribution (from 40 to 60 kHz per channel). Piecewise-linear GRAPE pulse is used. Calculation time: minutes

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

- Lines 19-20: 600 MHz magnet; implemented by `sys.magnet=14.1`.
- Lines 22-23: Isotopes; implemented by `sys.isotopes={'2H'}`.
- Lines 25-26: Alanine CD3 NQI parameters; implemented by `inter.coupling.matrix{1,1}=anas2mat(0,40e3,0,0,0,0)`.
- Lines 28-29: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 32-33: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 36-37: Powder context parameters; implemented by `parameters.spins={'2H'}`.
- Lines 44-45: Drift Liouvillians for the ensemble; implemented by `control.drifts=drifts(spin_system,@powder,parameters,'labframe')`.
- Lines 47-48: Set up and normalise the initial state; implemented by `rho_init=state(spin_system,'Lz','2H')`.
- Lines 51-52: Set up and normalise the target state; implemented by `rho_targ=state(spin_system,'Lx','2H')`.
- Lines 55-56: Get control operators; implemented by `Lx=operator(spin_system,'Lx','2H')`.
- Lines 59-60: Define control parameters; implemented by `control.isotopes={'2H'}`.
- Lines 75-76: Freeze the edges; implemented by `control.freeze=false(2,76)`.
- Lines 80-81: Plotting options; implemented by `control.plotting={'xy_controls','robustness','spectrogram'}`.
- Lines 83-84: Spinach housekeeping; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 86-87: Waveform guess with zeros at the edges; implemented by `guess=randn(2,76)/10; guess(:,1:4)=0`.
- Lines 90-91: Run the optimisation; implemented by `pulse_profile=fmaxnewton(spin_system,@grape_xy,guess)`.
- Lines 93-94: Get Cartesian components; implemented by `CLx=mean(control.pwr_levels)*pulse_profile(1,:)`.
- Lines 97-99: Apply RLC response distortion by a probe circuit with Q = 200; implemented by `kfigure(); restrans(CLx',CLy',control.pulse_dt(1), sys.magnet*spin('2H'),200,'pwl_tsc',100)`.

### Key state/data transformations

- Lines 20: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 23: computes `sys.isotopes` using `sys.isotopes={'2H'}`.
- Lines 26: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=anas2mat(0,40e3,0,0,0,0)`.
- Lines 29: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 30: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 33: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 37: computes `parameters.spins` using `parameters.spins={'2H'}`.
- Lines 38: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 39: computes `parameters.offset` using `parameters.offset=0`.
- Lines 40: computes `parameters.grid` using `parameters.grid='rep_2ang_100pts_sph'`.
- Lines 41: computes `parameters.rframes` using `parameters.rframes={{'2H',2}}`.
- Lines 42: computes `parameters.verbose` using `parameters.verbose=0`.
- Lines 45: computes `control.drifts` using `control.drifts=drifts(spin_system,@powder,parameters,'labframe')`.
- Lines 48: computes `rho_init` using `rho_init=state(spin_system,'Lz','2H')`.
- Lines 52: computes `rho_targ` using `rho_targ=state(spin_system,'Lx','2H')`.
- Lines 56: computes `Lx` using `Lx=operator(spin_system,'Lx','2H')`.
- Lines 57: computes `Ly` using `Ly=operator(spin_system,'Ly','2H')`.
- Lines 60: computes `control.isotopes` using `control.isotopes={'2H'}`.

## Implementation structure

- Probe circuit response effect on the accuracy of the deu-
- terium pre-phasing pulse designed to set deuterium magne-
- tisation in a -CD3 group of alanine up for rephasing 100
- microseconds after the pulse is finished.
- The system is assumed to be a powder (100 orientations)
- with a B1 distribution (from 40 to 60 kHz per channel).
- Piecewise-linear GRAPE pulse is used.
- Calculation time: minutes
- 600 MHz magnet
- Isotopes
- Alanine CD3 NQI parameters
- Basis set

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `anas2mat()`, `create()`, `basis()`, `drifts()`, `state()`, `operator()`, `false()`, `optimcon()`, `guess()`, `fmaxnewton()`, `pulse_profile()`, `kfigure()`, `restrans()`, `spin()`.
