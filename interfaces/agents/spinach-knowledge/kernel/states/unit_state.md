# kernel/states/unit_state.m

- Signature: `rho=unit_state(spin_system)`

## Purpose

Returns the unit state vector or matrix for the current formalism and basis. Syntax: `rho=unit_state(spin_system)`.

## Physical / mathematical content

The representation depends on the formalism: `sphten-liouv` returns the population of the `T(0,0)` basis state; `zeeman-liouv` returns the vectorised identity normalized to unit 2-norm; `zeeman-hilb` returns the identity matrix.

## Numerical / algorithmic content

The function selects a representation by `spin_system.bas.formalism`; unsupported values raise an error.

## Parameters / inputs

- `spin_system` - Spinach data object containing basis information; call `basis.m` first.

## Outputs

- `rho` - vector or matrix representation of the unit state.

## Implementation structure

- Checks that the basis formalism is present, then constructs the unit state in the corresponding representation.
