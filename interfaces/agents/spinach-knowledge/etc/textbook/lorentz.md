# etc/textbook/lorentz.m

- Signature: `[J,K,Kil]=lorentz(L)`

## Purpose

The (L,0)(+)(0,L) irreducible matrix representation of the Lorentz group with inversion. Syntax: [J,K,Kil]=lorentz(L)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- L -irreducible representation
- rank, e.g. 1/2

## Outputs

- J -three rotation generators
- K -three boost generators
- Kil -Killing form (expensive)

## Implementation structure

- The (L,0)(+)(0,L) irreducible matrix representation of the
- Lorentz group with inversion. Syntax:
- [J,K,Kil]=lorentz(L)
- L - irreducible representation
- rank, e.g. 1/2
- J - three rotation generators
- K - three boost generators
- Kil - Killing form (expensive)
- Check consistency
- Dimension and Pauli blocks
- Rotation and boost generators
- Collect generators
