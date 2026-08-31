# kernel/optimcon/tgrape.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/tgrape.m`
- Signature: `[fidelity,grad]=tgrape(spin_system,drift,controls,waveform,...`
- Total lines: 168

## Purpose

A special case of Gradient Ascent Pulse Engineering (GRAPE) objective function and gradient with respect to the vector of waveform slice du- rations. Syntax: [fidelity,grad]=tgrape(spin_system,drift,controls,waveform,... dt_grid,time_unit,rho_init,rho_targ)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 46-47: Check consistency; implemented by `grumble(spin_system,drift,controls,waveform,dt_grid,time_unit,rho_init,rho_targ)`.
- Lines 49-50: Count the time steps; implemented by `nsteps=size(waveform,2)`.
- Lines 52-53: Convert time units; implemented by `dt_grid=dt_grid*time_unit`.
- Lines 55-56: Preallocate trajectories; implemented by `fwd_traj=zeros([size(rho_init,1) (nsteps+1)],'like',1i)`.
- Lines 59-60: Hush up the output; implemented by `spin_system.sys.output='hush'`.
- Lines 62-63: Initialise forward and backward trajectories; implemented by `fwd_traj(:,1)=rho_init; bwd_traj(:,1)=rho_targ`.
- Lines 65-66: Compute trajectories; implemented by `for n=1:nsteps`.
- Lines 68-69: Start with the drift generator; implemented by `L_forw=drift; L_back=drift'`.
- Lines 71-72: Add current controls; implemented by `for k=1:numel(controls)`.
- Lines 77-78: Take time steps forwards and backwards; implemented by `fwd_traj(:,n+1)=step(spin_system,L_forw,fwd_traj(:,n),+dt_grid(n))`.
- Lines 83-84: Flip the backward trajectory; implemented by `bwd_traj=fliplr(bwd_traj)`.
- Lines 86-87: Compute Re(<a|b>) fidelity; implemented by `fidelity=real(rho_targ'*fwd_traj(:,end))`.
- Lines 89-90: Compute gradient; implemented by `if nargout>1`.
- Lines 92-93: Preallocate results; implemented by `grad=zeros(size(dt_grid))`.
- Lines 95-96: Loop over control sequence; implemented by `for n=1:nsteps`.
- Lines 98-99: Make evolution generator; implemented by `L=drift`.
- Lines 104-105: Compute fidelity derivative; implemented by `grad(n)=real(-1i*bwd_traj(:,n)'*L*fwd_traj(:,n))`.
- Lines 109-110: Convert time units; implemented by `grad=grad*time_unit`.

### Control flow inferred from the code

- Line 66: `for` loop over `n=1:nsteps`.
- Line 72: `for` loop over `k=1:numel(controls)`.
- Line 90: conditional branch on `nargout>1`.
- Line 96: `for` loop over `n=1:nsteps`.
- Line 100: `for` loop over `k=1:numel(controls)`.

### Key state/data transformations

- Lines 50: computes `nsteps` using `nsteps=size(waveform,2)`.
- Lines 53: computes `dt_grid` using `dt_grid=dt_grid*time_unit`.
- Lines 56: computes `fwd_traj` using `fwd_traj=zeros([size(rho_init,1) (nsteps+1)],'like',1i)`.
- Lines 57: computes `bwd_traj` using `bwd_traj=zeros([size(rho_init,1) (nsteps+1)],'like',1i)`.
- Lines 60: computes `spin_system.sys.output` using `spin_system.sys.output='hush'`.
- Lines 63: computes `fwd_traj(:,1)` using `fwd_traj(:,1)=rho_init; bwd_traj(:,1)=rho_targ`.
- Lines 69: computes `L_forw` using `L_forw=drift; L_back=drift'`.
- Lines 74: computes `L_back` using `L_back=L_back+waveform(k,nsteps+1-n)*controls{k}'`.
- Lines 78: computes `fwd_traj(:,n+1)` using `fwd_traj(:,n+1)=step(spin_system,L_forw,fwd_traj(:,n),+dt_grid(n))`.
- Lines 79: computes `bwd_traj(:,n+1)` using `bwd_traj(:,n+1)=step(spin_system,L_back,bwd_traj(:,n),-dt_grid(nsteps+1-n))`.
- Lines 87: computes `fidelity` using `fidelity=real(rho_targ'*fwd_traj(:,end))`.
- Lines 93: computes `grad` using `grad=zeros(size(dt_grid))`.
- Lines 99: computes `L` using `L=drift`.
- Lines 105: computes `grad(n)` using `grad(n)=real(-1i*bwd_traj(:,n)'*L*fwd_traj(:,n))`.

### Local helper functions

- Line 117: `grumble()` — `function grumble(spin_system,drift,controls,waveform,dt_grid,time_unit,rho_init,rho_targ)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})`.
  - Representative operation: `error('optimal control module requires Lioville space formalism.')`.

## Parameters / inputs

- drift -the drift Liouvillian, a matrix
- controls -control operators, a cell array
- of matrices
- waveform -control coefficients for each control
- operator (columns) at each time slice
- (rows), rad/s
- dt_grid -time slice durations, a col vector
- in the units of time chosen so that
- the elements are of the order of 1
- time_unit -unit of time, seconds; this is needed
- because optimisers get stuck when the
- variables are badly scaled
- rho_init -initial state of the system, a column
- vector
- rho_targ -target state of the system, a column
- vector

## Outputs

- fidelity -fidelity of the control sequence
- grad -gradient of the fidelity with respect to
- the durations of the slices

## Implementation structure

- A special case of Gradient Ascent Pulse Engineering (GRAPE) objective
- function and gradient with respect to the vector of waveform slice du-
- rations. Syntax:
- [fidelity,grad]=tgrape(spin_system,drift,controls,waveform,...
- dt_grid,time_unit,rho_init,rho_targ)
- drift -the drift Liouvillian, a matrix
- controls -control operators, a cell array
- of matrices
- waveform -control coefficients for each control
- operator (columns) at each time slice
- (rows), rad/s
- dt_grid -time slice durations, a col vector

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `fwd_traj()`, `bwd_traj()`, `waveform()`, `step()`, `dt_grid()`, `fliplr()`, `grad()`, `ismember()`, `iscolumn()`, `iscell()`, `any()`.
