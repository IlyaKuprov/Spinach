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

- Nonempty keyhole schedules with `newton` or `goodwin` are explicitly not implemented in `sphten-liouv` and `zeeman-liouv`, both in `optimcon` setup and direct `grape_liouv` calls. First-order keyhole methods, empty schedules, and existing Hilbert-space cases are unchanged; no algorithm is substituted.
- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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
