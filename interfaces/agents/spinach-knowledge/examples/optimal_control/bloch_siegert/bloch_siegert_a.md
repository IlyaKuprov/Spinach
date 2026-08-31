# examples/optimal_control/bloch_siegert/bloch_siegert_a.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/bloch_siegert/bloch_siegert_a.m`
- Signature: `bloch_siegert_a()`
- Total lines: 106

## Purpose

Bloch-Siegert shift compensation functionality demo. The script optimises a 90-degree pulse (Lz -> Lx) for a sing- le spin on resonance. As the control power is increased, Bloch-Siegert shift starts to reduce the fidelity unless it is correctly accounted for. Calculation time: minutes.

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

- Lines 13-14: Set magnetic field; implemented by `sys.magnet=14.1`.
- Lines 16-17: Set isotopes; implemented by `sys.isotopes={'1H'}`.
- Lines 19-20: Set interactions; implemented by `inter.zeeman.scalar={0.0}`.
- Lines 22-23: Set basis; implemented by `bas.formalism='sphten-liouv'`.
- Lines 26-27: Run Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 30-31: Build and normalise the initial state (Lz); implemented by `rho_init=state(spin_system,{'Lz'},{1})`.
- Lines 34-35: Build and normalise the target state (Lx); implemented by `rho_targ=state(spin_system,{'Lx'},{1})`.
- Lines 38-39: Get the control operators; implemented by `Lx=operator(spin_system,'Lx','1H')`.
- Lines 42-43: Build the drift Hamiltonian; implemented by `D=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 45-46: Optimal control settings; implemented by `control.isotopes={'1H'}`.
- Lines 57-58: Power levels to sweep (rad/s); implemented by `zeeman_frq=abs(spin('1H')*sys.magnet)`.
- Lines 61-62: Preallocate results; implemented by `fid_a=zeros(size(pwr_list))`.
- Lines 65-66: Common initial guess; implemented by `guess=randn(2,50)/10`.
- Lines 68-69: Loop over control power; implemented by `for n=1:numel(pwr_list)`.
- Lines 71-72: Set current power; implemented by `control.pwr_levels=pwr_list(n)`.
- Lines 74-75: Set pulse slice duration; implemented by `dt=(pi/100)/control.pwr_levels`.
- Lines 78-79: BSS settings: off, on; implemented by `control.bsiegert=false()`.
- Lines 84-85: Optimisation with and without BSS; implemented by `pulse_a=fmaxnewton(setting_a,@grape_xy,guess)`.

### Control flow inferred from the code

- Line 69: `for` loop over `n=1:numel(pwr_list)`.

### Key state/data transformations

- Lines 14: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 17: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 20: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.0}`.
- Lines 23: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 24: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 27: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 31: computes `rho_init` using `rho_init=state(spin_system,{'Lz'},{1})`.
- Lines 35: computes `rho_targ` using `rho_targ=state(spin_system,{'Lx'},{1})`.
- Lines 39: computes `Lx` using `Lx=operator(spin_system,'Lx','1H')`.
- Lines 40: computes `Ly` using `Ly=operator(spin_system,'Ly','1H')`.
- Lines 43: computes `D` using `D=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 46: computes `control.isotopes` using `control.isotopes={'1H'}`.
- Lines 47: computes `control.channels` using `control.channels=[1; 1]`.
- Lines 48: computes `control.drifts` using `control.drifts={{D}}`.
- Lines 49: computes `control.operators` using `control.operators={Lx,Ly}`.
- Lines 50: computes `control.rho_init` using `control.rho_init={rho_init}`.
- Lines 51: computes `control.rho_targ` using `control.rho_targ={rho_targ}`.
- Lines 52: computes `control.method` using `control.method='lbfgs'`.

## Implementation structure

- Bloch-Siegert shift compensation functionality demo. The
- script optimises a 90-degree pulse (Lz -> Lx) for a sing-
- le spin on resonance. As the control power is increased,
- Bloch-Siegert shift starts to reduce the fidelity unless
- it is correctly accounted for.
- Calculation time: minutes.
- Set magnetic field
- Set isotopes
- Set interactions
- Set basis
- Run Spinach housekeeping
- Build and normalise the initial state (Lz)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `hamiltonian()`, `assume()`, `spin()`, `pwr_list()`, `false()`, `optimcon()`, `true()`, `fmaxnewton()`, `fid_a()`, `ensemble()`, `fid_b()`, `kfigure()`.
