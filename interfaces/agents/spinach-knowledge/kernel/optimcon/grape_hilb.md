# kernel/optimcon/grape_hilb.m

- Signature: `[traj_data,fidelity,grad,hess]=grape_hilb(spin_system,drifts,controls,...`

## Purpose

Gradient Ascent Pulse Engineering (GRAPE) objective function, gradient and Hessian. Propagates the system through a user-supplied shaped pulse from a given initial state and projects the result onto the given final state. The fidelity is returned, along with its gradient and Hessian with respect to amplitudes of all control operators at every time step of the shaped pulse. Uses Hilbert-space formalism. Syntax: [traj_

## Physical / mathematical content

- Propagates Hilbert-space density matrices under drift and waveform-weighted control operators, and evaluates overlap-based fidelity with the target state.
- Supports piecewise-constant `rectangle` intervals and piecewise-linear `trapezium` intervals; the latter uses an Iserles-Nørsett product-quadrature generator.
- Fidelity can be evaluated at the terminal node or averaged over pulse nodes; configured trajectory-penalty operators are also averaged over those nodes.

## Numerical / algorithmic content

Zero fidelities and gradients are returned as valid values, including for auxiliary costates used by `grape_coop`. Initial-guess checks remain in `fmaxnewton`, where they apply to the assembled optimisation objective rather than individual GRAPE contributions.

- The routine precomputes interval propagators, applies them in a forward state sweep, and uses a backward costate sweep to form waveform derivatives.
- Hessians are available only with the `rectangle` integrator and `newton` or `goodwin` methods; trajectory cost terms do not support Hessians.

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
