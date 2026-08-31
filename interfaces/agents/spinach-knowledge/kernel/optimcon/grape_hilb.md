# kernel/optimcon/grape_hilb.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/grape_hilb.m`
- Signature: `[traj_data,fidelity,grad,hess]=grape_hilb(spin_system,drifts,controls,...`
- Total lines: 726

## Purpose

Gradient Ascent Pulse Engineering (GRAPE) objective function, gradient and Hessian. Propagates the system through a user-supplied shaped pulse from a given initial state and projects the result onto the given final state. The fidelity is returned, along with its gradient and Hessian with respect to amplitudes of all control operators at every time step of the shaped pulse. Uses Hilbert-space formalism. Syntax: [traj_

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 63-64: Check consistency; implemented by `grumble(spin_system,drifts,controls,waveform,rho_init,rho_targ,fidelity_type)`.
- Lines 66-67: Count the outputs; implemented by `n_outputs=nargout()`.
- Lines 69-70: Extract the timing grid; implemented by `dt=spin_system.control.pulse_dt`.
- Lines 72-73: Run array preallocations; implemented by `switch spin_system.control.integrator`.
- Lines 75-76: Piecewise-constant; implemented by `case 'rectangle'`.
- Lines 78-79: Number of time intervals and control operators; implemented by `nsteps=size(waveform,2); nctrls=size(waveform,1)`.
- Lines 81-82: Piecewise-linear; implemented by `case 'trapezium'`.
- Lines 84-85: Number of time intervals and control operators; implemented by `nsteps=size(waveform,2)-1; nctrls=size(waveform,1)`.
- Lines 89-90: Complain and bomb out; implemented by `error('unknown integrator: must be ''rectangle'' or ''trapezium''.')`.
- Lines 94-95: Hush up the output; implemented by `spin_system.sys.output='hush'`.
- Lines 97-99: Pull the target back through the dead time using the last drift generator in the drift array; implemented by `if spin_system.control.dead_time~=0`.
- Lines 104-106: Push the source through the prefix sequence using the first element of the drift array; implemented by `if ~isempty(spin_system.control.prefix)`.
- Lines 111-114: Push the target through the suffix sequence (which user needs to code in reverse time) using the last element of the drift array; implemented by `if ~isempty(spin_system.control.suffix)`.
- Lines 119-120: Preallocate forward trajectory; implemented by `fwd_traj=cell(1,nsteps+1); fwd_traj{1}=rho_init`.
- Lines 122-123: Preallocate backward trajectory; implemented by `if n_outputs>2`.
- Lines 127-128: Count the drifts; implemented by `ndrifts=numel(drifts)`.
- Lines 130-131: Pull Bloch-Siegert response operators; implemented by `bss_on=isfield(spin_system.control,'bsiegert')&&spin_system.control.bsiegert`.
- Lines 136-137: Precompute interval Hamiltonians and propagators; implemented by `switch spin_system.control.integrator`.

### Control flow inferred from the code

- Line 73: dispatches on `spin_system.control.integrator`; cases `'rectangle'`, `'trapezium'`.
- Line 99: conditional branch on `spin_system.control.dead_time~=0`.
- Line 106: conditional branch on `~isempty(spin_system.control.prefix)`.
- Line 114: conditional branch on `~isempty(spin_system.control.suffix)`.
- Line 123: conditional branch on `n_outputs>2`.
- Line 132: conditional branch on `bss_on`.
- Line 137: dispatches on `spin_system.control.integrator`; cases `'rectangle'`, `'trapezium'`.
- Line 146: `for` loop over `n=1:nsteps`.
- Line 152: `for` loop over `k=1:nctrls`.
- Line 160: conditional branch on `bss_on`.
- Line 161: `for` loop over `k=1:nctrls`.
- Line 178: `for` loop over `n=1:nsteps`.
- Line 185: `for` loop over `k=1:nctrls`.
- Line 208: `for` loop over `n=1:nsteps`.

### Key state/data transformations

- Lines 67: computes `n_outputs` using `n_outputs=nargout()`.
- Lines 70: computes `dt` using `dt=spin_system.control.pulse_dt`.
- Lines 79: computes `nsteps` using `nsteps=size(waveform,2); nctrls=size(waveform,1)`.
- Lines 95: computes `spin_system.sys.output` using `spin_system.sys.output='hush'`.
- Lines 100-101: computes `rho_targ` using `rho_targ=step(spin_system,drifts{end},rho_targ, -spin_system.control.dead_time)`.
- Lines 107: computes `prefix` using `prefix=spin_system.control.prefix`.
- Lines 108: computes `rho_init` using `rho_init=prefix(spin_system,drifts{1},rho_init)`.
- Lines 115: computes `suffix` using `suffix=spin_system.control.suffix`.
- Lines 120: computes `fwd_traj` using `fwd_traj=cell(1,nsteps+1); fwd_traj{1}=rho_init`.
- Lines 124: computes `bwd_traj` using `bwd_traj=cell(1,nsteps+1); bwd_traj{1}=rho_targ`.
- Lines 128: computes `ndrifts` using `ndrifts=numel(drifts)`.
- Lines 131: computes `bss_on` using `bss_on=isfield(spin_system.control,'bsiegert')&&spin_system.control.bsiegert`.
- Lines 133: computes `resp_ops` using `resp_ops=spin_system.control.resp_ops`.
- Lines 143: computes `H` using `H=cell(1,nsteps); P=cell(1,nsteps)`.
- Lines 149: computes `H{n}` using `H{n}=drifts{mod(n-1,ndrifts)+1}`.
- Lines 167: computes `P{n}` using `P{n}=propagator(spin_system,H{n},dt(n))`.
- Lines 181: computes `left_ham` using `left_ham=drifts{mod(n-1,ndrifts)+1}`.
- Lines 182: computes `right_ham` using `right_ham=drifts{mod(n,ndrifts)+1}`.

### Local helper functions

- Line 645: `grumble()` — `function grumble(spin_system,drifts,controls,waveform,rho_init,rho_targ,fidelity_type)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'zeeman-hilb'})`.
  - Representative operation: `error('this function requires a density matrix based formalism.')`.

## Parameters / inputs

- spin_system -Spinach data object that has been through
- the optimcon.m problem setup function.
- drifts -the drift Hamiltonians: a cell array con-
- taining one matrix (for time-independent
- drift) or multiple matrices (for time-de-
- pendent drift).
- controls -control operators in Hilbert space (cell
- array of matrices).
- waveform -control coefficients for each control ope-
- rator (in vertical dimension) at each time
- step (in horizonal dimension), rad/s
- rho_init -initial state of the system as a density
- matrix in Hilbert space.
- rho_targ -target state of the system as a density
- matrix in Hilbert space.
- fidelity_type -'real' (real part of the overlap)
- 'imag' (imaginary part of the overlap)
- 'square' (absolute square of the overlap)

## Outputs

- fidelity -fidelity of the control sequence
- grad -gradient of the fidelity with respect to
- the control sequence
- hess -Hessian of the fidelity with respect to the
- control sequence
- traj_data.forward -forward trajectory from the initial condi-
- tion(a stack of state matrices)
- Note: this is a low level function that is not designed to be called
- directly. Use grape_xy.m and grape_phase.m instead.

## Implementation structure

- Gradient Ascent Pulse Engineering (GRAPE) objective function, gradient
- and Hessian. Propagates the system through a user-supplied shaped pulse
- from a given initial state and projects the result onto the given final
- state. The fidelity is returned, along with its gradient and Hessian
- with respect to amplitudes of all control operators at every time step
- of the shaped pulse. Uses Hilbert-space formalism. Syntax:
- [traj_data,fidelity,grad,hess]=grape_hilb(spin_system,drifts,controls,...
- waveform,rho_init,rho_targ,...
- fidelity_type)
- spin_system -Spinach data object that has been through
- the optimcon.m problem setup function.
- drifts -the drift Hamiltonians: a cell array con-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `nargout()`, `step()`, `prefix()`, `suffix()`, `isfield()`, `waveform()`, `propagator()`, `keyhole_forw()`, `keyhole_back()`, `fliplr()`, `hdot()`, `dirdiff()`, `grad_col()`, `grad()`, `aux_mat()`.
