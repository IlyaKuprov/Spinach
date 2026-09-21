# kernel/operators/stevens.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/operators/stevens.m`
- Signature: `S=stevens(mult,k,q)`
- Total lines: 95

## Purpose

Extended Stevens operators. Syntax: S=stevens(mult,k,q)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- mult -multiplicity of the spin in question
- k -Stevens operator rank, a non-negative
- integer
- q -Stevens operator projection, an
- integer between -k and k

## Outputs

- S -Stevens operator matrix
- Note: for historical reasons, the definition of Stevens operators
- is irregular and must rely on explicitly stockpiled coeffi-
- cients. Only ranks smaller or equal to 12 are available.

## Implementation structure

- Extended Stevens operators. Syntax:
- S=stevens(mult,k,q)
- mult -multiplicity of the spin in question
- k -Stevens operator rank, a non-negative
- integer
- q -Stevens operator projection, an
- integer between -k and k
- S -Stevens operator matrix
- Note: for historical reasons, the definition of Stevens operators
- is irregular and must rely on explicitly stockpiled coeffi-
- cients. Only ranks smaller or equal to 12 are available.
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `pauli()`, `logical()`, `isscalar()`.
