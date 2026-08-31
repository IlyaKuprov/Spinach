# examples/optimal_control/state_transfer_coop.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/state_transfer_coop.m`
- Signature: `state_transfer_coop()`
- Total lines: 89

## Purpose

Optimal control pulse optimisation for state-to-state transfer in a quadrupolar 14N spin at a fixed orientation and power level Two pulses are optimised cooperatively, so that the sum of the outcomes only contains the target state, and no impurities. Calculation time: minutes.

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

- Lines 13-14: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 16-17: Isotopes; implemented by `sys.isotopes={'14N'}`.
- Lines 19-20: Glycine NQI, random orientation; implemented by `euler_angles=[1.0 2.0 3.0]`.
- Lines 23-24: Glycine 14N chemical shift; implemented by `inter.zeeman.scalar{1}=32.4`.
- Lines 26-27: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 30-31: Run Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 34-35: Set up and normalise the initial state; implemented by `rho_init=state(spin_system,'T1,0','14N')`.
- Lines 38-39: Set up and normalise the target state; implemented by `rho_targ=state(spin_system,'T2,0','14N')`.
- Lines 42-43: spin_system assumptions; implemented by `spin_system=assume(spin_system,'qnmr')`.
- Lines 45-46: Get the drift Hamiltonian; implemented by `[Iso,Q]=hamiltonian(spin_system)`.
- Lines 51-52: Get the control operators; implemented by `Lx=operator(spin_system,'Lx','14N')`.
- Lines 55-56: Define control parameters; implemented by `control.isotopes={'14N'}`.
- Lines 68-69: Plotting options; implemented by `control.plotting={'coherence_order','phi_controls'}`.
- Lines 71-72: Initial guess; implemented by `guess=randn(2,100)`.
- Lines 74-75: Spinach housekeeping; implemented by `spin_system=optimcon(spin_system,control)`.
- Lines 77-78: Optimisation; implemented by `[~,data]=fmaxnewton(spin_system,@grape_coop,guess)`.
- Lines 80-81: Pull out final states; implemented by `outcome_a=data.traj_data{1}{1}.forward(:,end)`.
- Lines 84-85: Print diagnostics; implemented by `disp(' Outcome A Outcome B (A+B)/2 Target')`.

### Key state/data transformations

- Lines 14: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 17: computes `sys.isotopes` using `sys.isotopes={'14N'}`.
- Lines 20: computes `euler_angles` using `euler_angles=[1.0 2.0 3.0]`.
- Lines 21: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=eeqq2nqi(1.18e6,0.53,1,euler_angles)`.
- Lines 24: computes `inter.zeeman.scalar{1}` using `inter.zeeman.scalar{1}=32.4`.
- Lines 27: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 28: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 31: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 35: computes `rho_init` using `rho_init=state(spin_system,'T1,0','14N')`.
- Lines 39: computes `rho_targ` using `rho_targ=state(spin_system,'T2,0','14N')`.
- Lines 46: computes `[Iso,Q]` using `[Iso,Q]=hamiltonian(spin_system)`.
- Lines 47: computes `H` using `H=Iso+orientation(Q,[1 2 3])`.
- Lines 48: computes `C` using `C=carrier(spin_system,'14N')`.
- Lines 52: computes `Lx` using `Lx=operator(spin_system,'Lx','14N')`.
- Lines 53: computes `Ly` using `Ly=operator(spin_system,'Ly','14N')`.
- Lines 56: computes `control.isotopes` using `control.isotopes={'14N'}`.
- Lines 57: computes `control.channels` using `control.channels=[1; 1]`.
- Lines 58: computes `control.drifts` using `control.drifts={{H}}`.

## Implementation structure

- Optimal control pulse optimisation for state-to-state transfer
- in a quadrupolar 14N spin at a fixed orientation and power level
- Two pulses are optimised cooperatively, so that the sum of the
- outcomes only contains the target state, and no impurities.
- Calculation time: minutes.
- Magnet field
- Isotopes
- Glycine NQI, random orientation
- Glycine 14N chemical shift
- Basis set
- Run Spinach housekeeping
- Set up and normalise the initial state

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `state()`, `assume()`, `hamiltonian()`, `orientation()`, `carrier()`, `rotframe()`, `operator()`, `optimcon()`, `fmaxnewton()`, `outcome_b()`.
