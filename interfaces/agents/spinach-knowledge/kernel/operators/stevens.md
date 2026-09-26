# kernel/operators/stevens.m

- Signature: `S=stevens(mult,k,q)`

## Purpose

Extended Stevens operators. Syntax: S=stevens(mult,k,q)

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

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
