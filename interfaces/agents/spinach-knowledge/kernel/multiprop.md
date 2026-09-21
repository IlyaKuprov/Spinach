# kernel/multiprop.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/multiprop.m`
- Signature: `rho=multiprop(spin_system,P,rho,N)`
- Total lines: 127

## Purpose

Applies a propagator repeatedly by binary adaptive squaring. Syntax: rho=multiprop(spin_system,P,rho,N)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- spin_system -Spinach spin system object
- P -propagator matrix
- rho -state vector or state-vector stack in Liouville space or
- wavefunction formalism, or a density matrix in Hilbert space
- formalism
- N -number of times to apply the propagator

## Outputs

- rho -state vector or density matrix after N applications of P
- Note: the algorithm expands N into binary powers, squares P successively,
- and applies only the active powers to rho. This avoids constructing
- P^N explicitly. Propagator squares are cleaned up using
- spin_system.tols.prop_chop.

## Implementation structure

- Applies a propagator repeatedly by binary adaptive squaring. Syntax:
- rho=multiprop(spin_system,P,rho,N)
- spin_system -Spinach spin system object
- P -propagator matrix
- rho -state vector or state-vector stack in Liouville space or
- wavefunction formalism, or a density matrix in Hilbert space
- formalism
- N -number of times to apply the propagator
- rho -state vector or density matrix after N applications of P
- Note: the algorithm expands N into binary powers, squares P successively,
- and applies only the active powers to rho. This avoids constructing
- P^N explicitly. Propagator squares are cleaned up using

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `uint64()`, `ismember()`, `bitand()`, `bitshift()`, `clean_up()`, `isstruct()`, `isfield()`, `ischar()`, `isscalar()`, `ismatrix()`, `isinteger()`, `allfinite()`.
