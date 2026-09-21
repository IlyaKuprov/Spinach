# kernel/operators/irr_sph_ten.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/operators/irr_sph_ten.m`
- Signature: `T=irr_sph_ten(mult,k)`
- Total lines: 108

## Purpose

Single-spin irreducible spherical tensor operators T(k,m) obeying the following commutation relation: [Lz,T_km]=m*T_km

## Physical / mathematical content

- Operator-construction utilities. They build bases and irreducible tensor representations for spin, bosonic, and transition operators.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Syntax

```matlab
T=irr_sph_ten(mult,k)
```

## Parameters / inputs

- mult -multiplicity of the spin in question
- k -irreducible spherical tensor rank (optional)

## Outputs

- T -a two-argument call returns a cell array of
- tensors of rank k in the order of decreasing
- projection. A single argument call produces
- tensors of all ranks and puts them into a
- cell array in the order of increasing rank,
- and decreasing projection within each rank.
- Note: operator normalisation in spin dynamics is not a good
- idea. The only way to make the formalism independent of the
- spin quantum number is to impose identical commutation rela-
- tions rather than equal matrix norms.

## Implementation structure

- Single-spin irreducible spherical tensor operators T(k,m)
- obeying the following commutation relation:
- [Lz,T_km]=m*T_km
- T=irr_sph_ten(mult,k)
- mult -multiplicity of the spin in question
- k -irreducible spherical tensor rank (optional)
- T -a two-argument call returns a cell array of
- tensors of rank k in the order of decreasing
- projection. A single argument call produces
- tensors of all ranks and puts them into a
- cell array in the order of increasing rank,
- and decreasing projection within each rank.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `speye()`, `pauli()`, `isscalar()`.
