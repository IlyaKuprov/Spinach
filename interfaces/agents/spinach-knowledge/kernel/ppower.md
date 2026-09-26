# kernel/ppower.m

- Signature: `P=ppower(spin_system,P,N)`

## Purpose

Computes integer propagator powers via an efficient powers-of-two based strategy. Syntax: P=ppower(spin_system,P,N)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- spin_system -Spinach spin system object
- P -propagator matrix
- N -non-negative integer propagator power

## Outputs

- P -propagator matrix raised to the power of N
- Note: the algorithm expands N into binary powers, squares P succes-
- sively, and multiplies only the active powers into the result.
- This avoids explicit repeated multiplication. Propagator pow-
- ers are cleaned up using spin_system.tols.prop_chop.

## Implementation structure

- Computes integer propagator powers via an efficient powers-of-two
- based strategy. Syntax:
- P=ppower(spin_system,P,N)
- spin_system -Spinach spin system object
- P -propagator matrix
- N -non-negative integer propagator power
- P -propagator matrix raised to the power of N
- Note: the algorithm expands N into binary powers, squares P succes-
- sively, and multiplies only the active powers into the result.
- This avoids explicit repeated multiplication. Propagator pow-
- ers are cleaned up using spin_system.tols.prop_chop.
- Check consistency
