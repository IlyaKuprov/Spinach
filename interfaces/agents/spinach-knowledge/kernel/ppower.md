# kernel/ppower.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/ppower.m`
- Signature: `P=ppower(spin_system,P,N)`
- Total lines: 113

## Purpose

Computes integer propagator powers via an efficient powers-of-two based strategy. Syntax: P=ppower(spin_system,P,N)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `uint64()`, `issparse()`, `speye()`, `bitand()`, `clean_up()`, `bitshift()`, `isstruct()`, `isfield()`, `isscalar()`, `ismatrix()`, `isinteger()`, `allfinite()`.
