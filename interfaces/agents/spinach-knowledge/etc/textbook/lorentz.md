# etc/textbook/lorentz.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/textbook/lorentz.m`
- Signature: `[J,K,Kil]=lorentz(L)`
- Total lines: 83

## Purpose

The (L,0)(+)(0,L) irreducible matrix representation of the Lorentz group with inversion. Syntax: [J,K,Kil]=lorentz(L)

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `pauli()`, `transpose()`, `Kil()`, `isscalar()`.
