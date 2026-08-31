# examples/optimal_control/features_bss.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/features_bss.m`
- Signature: `features_bss()`
- Total lines: 90

## Purpose

Optimal control pulse optimisation with Bloch-Siegert shift corrections switched on. A single proton with a Larmor frequency of 1 MHz is driven at a significant fraction of its Larmor frequency --a regime where the counter-rotating component of the control field shifts the resonance ap- preciably. A 90-degree pulse is optimised using LBFGS-GRAPE algorithm. For comparison, the same pulse is also optimised with the cor

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

- Lines 17-18: Larmor frequency of 1 MHz; implemented by `sys.magnet=2*pi*1e6/spin('1H')`.
- Lines 20-21: Spin system; implemented by `sys.isotopes={'1H'}`.
- Lines 23-24: Chemical shifts, ppm; implemented by `inter.zeeman.scalar={0.0}`.
- Lines 26-27: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 30-31: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 34-35: Set up and normalise the initial state; implemented by `rho_init=state(spin_system,{'Lz'},{1})`.
- Lines 38-39: Set up and normalise the target state; implemented by `rho_targ=state(spin_system,{'Lx'},{1})`.
- Lines 42-43: Drift Liouvillian is zero on resonance; implemented by `dim=size(spin_system.bas.basis,1)`.
- Lines 45-46: Control parameters; implemented by `control.drifts={{sparse(dim,dim)}}`.
- Lines 61-62: Deterministic initial guess; implemented by `guess=[linspace(0.1,0.5,50); 0.05*ones(1,50)]`.
- Lines 64-65: Optimisation; implemented by `control.bsiegert=true()`.
- Lines 69-70: Corrected pulse in the corrected model; implemented by `[~,fid_corr]=ensemble(pulse,spin_system)`.
- Lines 72-73: Optimisation; implemented by `control.bsiegert=false()`.
- Lines 77-78: Uncorrected pulse in the corrected model; implemented by `control.bsiegert=true()`.
- Lines 82-83: Report the results; implemented by `report(spin_system,'---------------')`.

### Key state/data transformations

- Lines 18: computes `sys.magnet` using `sys.magnet=2*pi*1e6/spin('1H')`.
- Lines 21: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 24: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.0}`.
- Lines 27: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 28: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 31: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 35: computes `rho_init` using `rho_init=state(spin_system,{'Lz'},{1})`.
- Lines 39: computes `rho_targ` using `rho_targ=state(spin_system,{'Lx'},{1})`.
- Lines 43: computes `dim` using `dim=size(spin_system.bas.basis,1)`.
- Lines 46: computes `control.drifts` using `control.drifts={{sparse(dim,dim)}}`.
- Lines 47-48: computes `control.operators` using `control.operators={operator(spin_system,'Lx','1H'), operator(spin_system,'Ly','1H')}`.
- Lines 49: computes `control.off_ops` using `control.off_ops={operator(spin_system,'Lx','1H')}`.
- Lines 50: computes `control.offsets` using `control.offsets={linspace(-1e5,1e5,11)}`.
- Lines 51: computes `control.isotopes` using `control.isotopes={'1H'}`.
- Lines 52: computes `control.channels` using `control.channels=[1;1]`.
- Lines 53: computes `control.rho_init` using `control.rho_init={rho_init}`.
- Lines 54: computes `control.rho_targ` using `control.rho_targ={rho_targ}`.
- Lines 55: computes `control.pwr_levels` using `control.pwr_levels=0.2*abs(spin_system.inter.basefrqs(1))`.

## Implementation structure

- Optimal control pulse optimisation with Bloch-Siegert shift corrections
- switched on. A single proton with a Larmor frequency of 1 MHz is driven
- at a significant fraction of its Larmor frequency --a regime where the
- counter-rotating component of the control field shifts the resonance ap-
- preciably. A 90-degree pulse is optimised using LBFGS-GRAPE algorithm.
- For comparison, the same pulse is also optimised with the corrections
- switched off and then evaluated in the corrected model.
- Calculation time: minutes.
- Larmor frequency of 1 MHz
- Spin system
- Chemical shifts, ppm
- Basis set

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `spin()`, `create()`, `basis()`, `state()`, `operator()`, `true()`, `optimcon()`, `fmaxnewton()`, `ensemble()`, `false()`, `report()`, `num2str()`.
