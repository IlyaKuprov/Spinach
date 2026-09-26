# kernel/states/unit_state.m

- Signature: `rho=unit_state(spin_system)`

## Purpose

Returns the unit state vector or matrix in the current formalism and basis. Syntax: rho=unit_state(spin_system)

## Physical / mathematical content

- State-construction utilities. These routines build equilibrium states, singlets, triplets, partner-state expansions, and physically meaningful density operators in the active basis.

## Numerical / algorithmic content

## Parameters / inputs

- spin_system -Spinach data object containing basis
- information (call basis.m first)

## Outputs

- rho -vector or matrix representation of
- the unit state

## Implementation structure

- Returns the unit state vector or matrix in the current formalism
- and basis. Syntax:
- rho=unit_state(spin_system)
- spin_system -Spinach data object containing basis
- information (call basis.m first)
- rho -vector or matrix representation of
- the unit state
- Check consistency
- Decide how to proceed
- Unit population of T(0,0) state
- Normalized stretched unit matrix
- Sparse unit matrix
