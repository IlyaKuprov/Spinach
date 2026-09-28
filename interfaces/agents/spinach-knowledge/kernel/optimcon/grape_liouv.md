# kernel/optimcon/grape_liouv.m

- Signature: `[traj_data,fidelity,grad,hess]=grape_liouv(spin_system,drifts,controls,...`

## Purpose

Gradient Ascent Pulse Engineering (GRAPE) objective function, gradient and Hessian. Propagates the system through a user-supplied shaped pulse from a given initial state and projects the result onto the given final state. The fidelity is returned, along with its gradient and Hessian with respect to amplitudes of all control operators at every time step of the shaped pulse. Uses Liouville-space or wavefunction formalisms.

## Numerical / algorithmic content

- The fidelity and its derivatives are propagated through the pulse sequence; both rectangular (piecewise-constant) and trapezium (piecewise-linear) integrators are supported.
- Nonempty keyhole schedules with `newton` or `goodwin` are explicitly not implemented in `sphten-liouv`, `zeeman-liouv`, and `zeeman-wavef`, both in `optimcon` setup and direct `grape_liouv` calls. First-order `lbfgs`/`rbfgs` keyhole methods and empty schedules remain available; no algorithm is substituted. A four-output Hessian request with keyholes is refused.
- Zero fidelities and gradients are valid, including for auxiliary costates used by `grape_coop`. Initial-guess checks belong in `fmaxnewton`, where they apply to the assembled optimisation objective rather than individual GRAPE contributions. Hessians are unavailable for piecewise-linear integration, stroboscopic steady states, nonempty keyhole schedules, or trajectory cost terms.

## Parameters / inputs

- spin_system -Spinach data object that has been through
- the optimcon.m problem setup function.
- drifts -drift generators (Liouvillians or wavefunction Hamiltonians):
- a cell array containing one matrix (for time-independent
- drift) or multiple matrices (one per time
- slice / point, for time-dependent drift).
- controls -control generators in the selected formalism (cell array of matrices).
- waveform -control coefficients for each control ope-
- rator (in vertical dimension) at each time
- slice / point (horizonal dimension), rad/s
- rho_init -initial state as a Liouville-space vector or wavefunction,
- ignored in stroboscopic
- steady state optimisations
- rho_targ -target state as a Liouville-space vector or wavefunction.
- fidelity_type -'real' (real part of the overlap)
- 'imag' (imaginary part of the overlap)
- 'square' (absolute square of the overlap)

## Outputs

- fidelity -fidelity of the control sequence
- grad -gradient of the fidelity with respect to
- the control sequence
- hess -Hessian of the fidelity with respect to
- the control sequence, not available for
- piecewise-linear or stroboscopic steady-state optimisations,
- or nonempty keyhole schedules
- traj_data.forward -forward trajectory from the initial con-
- dition or stroboscopic steady state (a
- stack of state vectors)
- Note: this is a low level function that is not designed to be called
- directly. Use grape_xy.m, grape_phase.m, or other wrapper func-
- tions instead.
- TODO (Keitel): add logic to avoid computing backward trajectory
- when the gradient is not requested
