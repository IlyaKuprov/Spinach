# kernel/optimcon/grape_hilb.m

- Signature: `[traj_data,fidelity,grad,hess]=grape_hilb(spin_system,drifts,controls,...`

## Purpose

Gradient Ascent Pulse Engineering (GRAPE) objective function, gradient and Hessian. Propagates the system through a user-supplied shaped pulse from a given initial state and projects the result onto the given final state. The fidelity is returned, along with its gradient and Hessian with respect to amplitudes of all control operators at every time step of the shaped pulse. Uses Hilbert-space formalism. Syntax: [traj_

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

Zero fidelities and gradients are returned as valid values, including for auxiliary costates used by `grape_coop`. Initial-guess checks remain in `fmaxnewton`, where they apply to the assembled optimisation objective rather than individual GRAPE contributions.

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

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
