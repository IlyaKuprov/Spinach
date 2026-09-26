# kernel/optimcon/grape_liouv.m

- Signature: `[traj_data,fidelity,grad,hess]=grape_liouv(spin_system,drifts,controls,...`

## Purpose

Gradient Ascent Pulse Engineering (GRAPE) objective function, gradient and Hessian. Propagates the system through a user-supplied shaped pulse from a given initial state and projects the result onto the given final state. The fidelity is returned, along with its gradient and Hessian with respect to amplitudes of all control operators at every time step of the shaped pulse. Uses Liouville-space or wavefunction formalisms.

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- Nonempty keyhole schedules with `newton` or `goodwin` are explicitly not implemented in `sphten-liouv`, `zeeman-liouv`, and `zeeman-wavef`, both in `optimcon` setup and direct `grape_liouv` calls. First-order `lbfgs`/`rbfgs` keyhole methods, empty schedules, and existing Hilbert-space keyhole Hessians remain available; no algorithm is substituted. The method restriction applies regardless of output count, and any four-output request with a state-vector keyhole is refused.
- Zero fidelities and gradients are returned as valid values, including for auxiliary costates used by `grape_coop`. Initial-guess checks remain in `fmaxnewton`, where they apply to the assembled optimisation objective rather than individual GRAPE contributions.

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

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
