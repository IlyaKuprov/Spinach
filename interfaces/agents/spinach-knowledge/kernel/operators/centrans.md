# kernel/operators/centrans.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/operators/centrans.m`
- Signature: `A=centrans(mult,type)`
- Total lines: 87

## Purpose

Central transition operators of half-integer spins in the Pauli basis. Syntax: A=centrans(mult,type)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- mult -multipicity of the spin in question, an
- even positive integer
- type -operator type: 'z' for polarisation, '+'
- for raising, '-' for lowering

## Outputs

- A -central transition operator

## Implementation structure

- Central transition operators of half-integer spins in the
- Pauli basis. Syntax:
- A=centrans(mult,type)
- mult -multipicity of the spin in question, an
- even positive integer
- type -operator type: 'z' for polarisation, '+'
- for raising, '-' for lowering
- A -central transition operator
- Check consistency
- Build CT operator
- Sx on central transition
- Sy on central transition

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spalloc()`, `complex()`, `isscalar()`, `ischar()`, `ismember()`.
