# examples/optimal_control/case_studies/Tosner_JMR_2009/coherence_transfer.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/case_studies/Tosner_JMR_2009/coherence_transfer.m`
- Signature: `coherence_transfer()`
- Total lines: 79

## Purpose

The first optimal control example from A heteronuclear two-spin system (1H–13C) with an scalar and both nuclei set on resonance; the goal is to trans- fer transverse magnetisation from proton to carbon: Hx → Cx over a fixed evolution period T = 1/J.

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

- Lines 17-18: Magnetic field, Tesla; implemented by `sys.magnet=14.1`.
- Lines 21-22: Chemical shifts, ppm; implemented by `inter.zeeman.scalar={0.0,0.0}`.
- Lines 24-25: Scalar coupling, Hz; implemented by `inter.coupling.scalar=cell(2,2)`.
- Lines 28-29: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 32-33: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 36-37: Initial state: Lx on proton (spin 1); implemented by `rho_init=state(spin_system,{'Lx'},{1})`.
- Lines 40-41: Target state: Lx on carbon (spin 2); implemented by `rho_targ=state(spin_system,{'Lx'},{2})`.
- Lines 44-45: Control operators; implemented by `LxH=operator(spin_system,'Lx',1)`.
- Lines 50-51: Drift Hamiltonian; implemented by `H=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 53-54: Control data structure; implemented by `control.isotopes={'1H','13C'}`.
- Lines 67-68: Visual diagnostics; implemented by `control.plotting={'xy_controls','spectrogram','robustness'}`.
- Lines 70-71: Random initial guess; implemented by `guess=rand(4,150)/10`.
- Lines 73-74: Optimisation; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 77-79: Better plotting needed here; implemented by `end`.

### Key state/data transformations

- Lines 18: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 19: computes `sys.isotopes` using `sys.isotopes={'1H','13C'}`.
- Lines 22: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.0,0.0}`.
- Lines 25: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(2,2)`.
- Lines 26: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=140`.
- Lines 29: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 30: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 33: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 37: computes `rho_init` using `rho_init=state(spin_system,{'Lx'},{1})`.
- Lines 41: computes `rho_targ` using `rho_targ=state(spin_system,{'Lx'},{2})`.
- Lines 45: computes `LxH` using `LxH=operator(spin_system,'Lx',1)`.
- Lines 46: computes `LyH` using `LyH=operator(spin_system,'Ly',1)`.
- Lines 47: computes `LxC` using `LxC=operator(spin_system,'Lx',2)`.
- Lines 48: computes `LyC` using `LyC=operator(spin_system,'Ly',2)`.
- Lines 51: computes `H` using `H=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 54: computes `control.isotopes` using `control.isotopes={'1H','13C'}`.
- Lines 55: computes `control.channels` using `control.channels=[1; 1; 2; 2]`.
- Lines 56: computes `control.drifts` using `control.drifts={{H}}`.

## Implementation structure

- The first optimal control example from
- A heteronuclear two-spin system (1H–13C) with an scalar
- and both nuclei set on resonance; the goal is to trans-
- fer transverse magnetisation from proton to carbon:
- Hx → Cx
- over a fixed evolution period T = 1/J.
- Magnetic field, Tesla
- Chemical shifts, ppm
- Scalar coupling, Hz
- Basis set
- Spinach housekeeping
- Initial state: Lx on proton (spin 1)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `hamiltonian()`, `assume()`, `optimcon()`, `fmaxnewton()`.
