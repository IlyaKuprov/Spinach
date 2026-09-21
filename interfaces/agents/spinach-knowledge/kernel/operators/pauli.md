# kernel/operators/pauli.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/operators/pauli.m`
- Signature: `S=pauli(mult)`
- Total lines: 98

## Purpose

Pauli spin operators (sparse, see below for normalisa- tion conventions) for a spin of a user-specified ener- gy level multiplicity. Syntax: S=pauli(mult)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- mult -an integer specifying the
- multiplicity of the spin

## Outputs

- S.u -unit operator
- S.p -raising operator
- S.m -lowering operator
- S.x -Sx observable operator
- S.y -Sy observable operator
- S.z -Sz observable operator
- Note: the matrices are normalised to obey the following
- commutation relations for all multiplicities:
- [S.x,S.y]=1i*S.z
- [S.y,S.z]=1i*S.x
- [S.z,S.x]=1i*S.y
- Note: raising and lowering operators are defined as:
- S.p=S.x+1i*S.y
- S.m=S.x-1i*S.y
- Note: arrays are declared complex at creation to avoid
- expensive reallocation operations later on.

## Implementation structure

- Pauli spin operators (sparse, see below for normalisa-
- tion conventions) for a spin of a user-specified ener-
- gy level multiplicity. Syntax:
- S=pauli(mult)
- mult -an integer specifying the
- multiplicity of the spin
- S.u -unit operator
- S.p -raising operator
- S.m -lowering operator
- S.x -Sx observable operator
- S.y -Sy observable operator
- S.z -Sz observable operator

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `double()`, `complex()`, `spdiags()`, `speye()`, `isscalar()`.
