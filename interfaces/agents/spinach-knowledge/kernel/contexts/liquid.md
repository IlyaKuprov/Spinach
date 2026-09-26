# kernel/contexts/liquid.m

- Signature: `answer=liquid(spin_system,pulse_sequence,parameters,assumptions)`

## Purpose

Liquid-phase interface to pulse sequences. Generates a Liouvillian superoperator and passes it on to the pulse sequence function, which should be supplied as a handle. This interface handles RDC mode --if the 'rdc' need is specified, it would use the order matrix supplied by the user to compute the residual anisotropies of all interactions. Syntax: answer=liquid(spin_system,pulse_sequence,parameters,assumptions)

## Physical / mathematical content

- Simulation-context constructors. These wrappers assemble Hamiltonians, Liouvillians, relaxation, kinetics, quadrature grids, and orientation/spatial machinery for a particular physical regime.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Parameters / inputs

- pulse_sequence -pulse sequence function handle. See the
- experiments directory for the list of
- pulse sequences that ship with Spinach.
- parameters.spins -a cell array giving the
- spins that the pulse sequence works on, in
- the order of channels, e.g. {'1H','13C'}
- parameters.offset -a cell array giving
- transmitter offsets on each of the spins
- listed in parameters.spins array.
- parameters.needs -a cell array of strings specifying additional
- information required by the sequence:
- 'zeeman_op' -Zeeman part of the Hamiltonian
- in the laboratory frame, to be placed into
- parameters.hzeeman and sent to pulse sequence
- 'rdc' -triggers the processing of residual
- anisotropic couplings due to partial order
- 'rho_eq' -thermal equilibrium state at the
- specified temperature with respect to the
- isotropic part of the Hamiltonian; this is
- placed into parameters.rho0
- parameters.rframes -rotating frame specification, e.g.
- {{'13C',2},{'14N,3}} requests second
- order rotating frame transformation
- with respect to carbon-13 and third
- order rotating frame transformation
- with respect to nitrogen-14. When
- this option is used, the assumptions
- on the respective spins should be
- laboratory frame.
- parameters.* -additional subfields may be required by
- the pulse sequence -check its docs
- assumptions -context-specific assumptions ('nmr', 'epr',
- 'labframe', etc.) -see the pulse sequence
- header for information on this setting.

## Outputs

- answer -whatever it is that the pulse sequence returns.
- Note: arbitrary order rotating frame transformation is supported, inc-
- luding infinite order. See the header of rotframe.m for further
- information.

## Implementation structure

- Liquid-phase interface to pulse sequences. Generates a Liouvillian
- superoperator and passes it on to the pulse sequence function, which
- should be supplied as a handle.
- This interface handles RDC mode --if the 'rdc' need is specified,
- it would use the order matrix supplied by the user to compute the
- residual anisotropies of all interactions. Syntax:
- answer=liquid(spin_system,pulse_sequence,parameters,assumptions)
- pulse_sequence -pulse sequence function handle. See the
- experiments directory for the list of
- pulse sequences that ship with Spinach.
- parameters.spins -a cell array giving the
- spins that the pulse sequence works on, in
