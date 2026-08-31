# kernel/optimcon/grape_liouv.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/grape_liouv.m`
- Signature: `[traj_data,fidelity,grad,hess]=grape_liouv(spin_system,drifts,controls,...`
- Total lines: 986

## Purpose

Gradient Ascent Pulse Engineering (GRAPE) objective function, gradient and Hessian. Propagates the system through a user-supplied shaped pulse from a given initial state and projects the result onto the given final state. The fidelity is returned, along with its gradient and Hessian with respect to amplitudes of all control operators at every time step of the shaped pulse. Uses Liouville-space formalism. Syntax: [tra

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 73-75: Check consistency; implemented by `grumble(spin_system,drifts,controls,waveform, rho_init,rho_targ,fidelity_type)`.
- Lines 77-78: Count the outputs; implemented by `n_outputs=nargout()`.
- Lines 80-81: Extract the timing grid; implemented by `dt=spin_system.control.pulse_dt`.
- Lines 83-84: Make freeze mask explicit; implemented by `if isempty(spin_system.control.freeze)`.
- Lines 90-91: Run array preallocations; implemented by `switch spin_system.control.integrator`.
- Lines 93-94: Piecewise-constant; implemented by `case 'rectangle'`.
- Lines 96-97: Number of time intervals and control operators; implemented by `nsteps=size(waveform,2); nctrls=size(waveform,1)`.
- Lines 99-100: Preallocate forward and backward trajectories; implemented by `fwd_traj=zeros([size(rho_init,1) (nsteps+1)],'like',1i)`.
- Lines 103-104: Preallocate arrays used in Hessian calculation; implemented by `if n_outputs>3`.
- Lines 106-107: Both Newton and Goodwin; implemented by `fwd_dP=cell(nctrls,nsteps)`.
- Lines 111-113: Goodwin Hessian route needs cumulative propagators; implemented by `if strcmp(spin_system.control.method,'goodwin')`.
- Lines 119-120: Empty declarations tested downstream; implemented by `fwd_dP={}; bwd_dP={}; fwd_d2P={}; P_cum={}`.
- Lines 124-125: Steady state needs propagators; implemented by `if spin_system.control.steady`.
- Lines 131-132: Piecewise-linear; implemented by `case 'trapezium'`.
- Lines 134-135: Number of time intervals and control operators; implemented by `nsteps=size(waveform,2)-1; nctrls=size(waveform,1)`.
- Lines 141-142: Cannot do trapezium Hessians yet; implemented by `fwd_dP={}; bwd_dP={}; fwd_d2P={}; P_cum={}`.
- Lines 146-147: Complain and bomb out; implemented by `error('unknown integrator: must be ''rectangle'' or ''trapezium''.')`.
- Lines 151-152: Hush up the output; implemented by `spin_system.sys.output='hush'`.

### Control flow inferred from the code

- Line 84: conditional branch on `isempty(spin_system.control.freeze)`.
- Line 91: dispatches on `spin_system.control.integrator`; cases `'rectangle'`, `'trapezium'`.
- Line 104: conditional branch on `n_outputs>3`.
- Line 113: conditional branch on `strcmp(spin_system.control.method,'goodwin')`.
- Line 125: conditional branch on `spin_system.control.steady`.
- Line 156: conditional branch on `spin_system.control.dead_time~=0`.
- Line 163: conditional branch on `~isempty(spin_system.control.prefix)`.
- Line 171: conditional branch on `~isempty(spin_system.control.suffix)`.
- Line 188: conditional branch on `bss_on`.
- Line 193: dispatches on `spin_system.control.integrator`; cases `'rectangle'`, `'trapezium'`.
- Line 202: `for` loop over `n=1:nsteps`.
- Line 211: `for` loop over `k=1:nctrls`.
- Line 222: conditional branch on `bss_on`.
- Line 223: `for` loop over `k=1:nctrls`.

### Key state/data transformations

- Lines 78: computes `n_outputs` using `n_outputs=nargout()`.
- Lines 81: computes `dt` using `dt=spin_system.control.pulse_dt`.
- Lines 85: computes `frozen` using `frozen=false(size(waveform))`.
- Lines 97: computes `nsteps` using `nsteps=size(waveform,2); nctrls=size(waveform,1)`.
- Lines 100: computes `fwd_traj` using `fwd_traj=zeros([size(rho_init,1) (nsteps+1)],'like',1i)`.
- Lines 101: computes `bwd_traj` using `bwd_traj=zeros([size(rho_init,1) (nsteps+1)],'like',1i)`.
- Lines 107: computes `fwd_dP` using `fwd_dP=cell(nctrls,nsteps)`.
- Lines 108: computes `bwd_dP` using `bwd_dP=cell(nctrls,nsteps)`.
- Lines 109: computes `fwd_d2P` using `fwd_d2P=cell(nctrls,nctrls,nsteps)`.
- Lines 114: computes `P_cum` using `P_cum=cell(1,nsteps)`.
- Lines 126: computes `P` using `P=cell(1,nsteps)`.
- Lines 152: computes `spin_system.sys.output` using `spin_system.sys.output='hush'`.
- Lines 157-158: computes `rho_targ` using `rho_targ=step(spin_system,drifts{end}',rho_targ, -spin_system.control.dead_time)`.
- Lines 164: computes `prefix` using `prefix=spin_system.control.prefix`.
- Lines 165: computes `rho_init` using `rho_init=prefix(spin_system,drifts{1},rho_init)`.
- Lines 172: computes `suffix` using `suffix=spin_system.control.suffix`.
- Lines 177: computes `zero_state` using `zero_state=complex(spalloc(size(rho_init,1),size(rho_init,2),0))`.
- Lines 178: computes `zero_drift` using `zero_drift=complex(spalloc(size(drifts{1},1),size(drifts{1},2),0))`.

### Local helper functions

- Line 905: `grumble()` — `function grumble(spin_system,drifts,controls,waveform,`.
  - Representative operation: `rho_init,rho_targ,fidelity_type)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv', 'zeeman-liouv', 'zeeman-wavef'})`.

## Parameters / inputs

- spin_system -Spinach data object that has been through
- the optimcon.m problem setup function.
- drifts -the drift Liouvillians: a cell array con-
- taining one matrix (for time-independent
- drift) or multiple matrices (one per time
- slice / point, for time-dependent drift).
- controls -control operators in Liouville space (cell
- array of matrices).
- waveform -control coefficients for each control ope-
- rator (in vertical dimension) at each time
- slice / point (horizonal dimension), rad/s
- rho_init -initial state of the system as a vector in
- Liouville space, ignored in stroboscopic
- steady state optimisations
- rho_targ -target state of the system as a vector in
- Liouville space.
- fidelity_type -'real' (real part of the overlap)
- 'imag' (imaginary part of the overlap)
- 'square' (absolute square of the overlap)

## Outputs

- fidelity -fidelity of the control sequence
- grad -gradient of the fidelity with respect to
- the control sequence
- hess -Hessian of the fidelity with respect to
- the control sequence, not available for
- piecewise-linear and stroboscopic stea-
- dy state optimisations
- traj_data.forward -forward trajectory from the initial con-
- dition or stroboscopic steady state (a
- stack of state vectors)
- Note: this is a low level function that is not designed to be called
- directly. Use grape_xy.m, grape_phase.m, or other wrapper func-
- tions instead.
- TODO (Keitel): add logic to avoid computing backward trajectory
- when the gradient is not requested

## Implementation structure

- Gradient Ascent Pulse Engineering (GRAPE) objective function, gradient
- and Hessian. Propagates the system through a user-supplied shaped pulse
- from a given initial state and projects the result onto the given final
- state. The fidelity is returned, along with its gradient and Hessian
- with respect to amplitudes of all control operators at every time step
- of the shaped pulse. Uses Liouville-space formalism. Syntax:
- [traj_data,fidelity,...
- grad,hess]=grape_liouv(spin_system,drifts,controls,...
- waveform,rho_init,rho_targ,...
- fidelity_type)
- spin_system -Spinach data object that has been through
- the optimcon.m problem setup function.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `nargout()`, `false()`, `strcmp()`, `step()`, `prefix()`, `suffix()`, `complex()`, `spalloc()`, `fwd_traj()`, `bwd_traj()`, `isfield()`, `waveform()`, `speye()`, `propagator()`, `clean_up()`.
