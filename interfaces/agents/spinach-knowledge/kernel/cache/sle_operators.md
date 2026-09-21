# kernel/cache/sle_operators.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/cache/sle_operators.m`
- Signature: `[Lx,Ly,Lz,D,space_basis]=sle_operators(max_rank,int_ranks)`
- Total lines: 347

## Purpose

Wigner D function basis set and rotation generators required by the SLE module. Syntax: [Lx,Ly,Lz,D,space_basis]=sle_operators(max_rank,int_ranks)

## Physical / mathematical content

- Cache-management utilities. These files maintain Spinach temporary or persistent cache state used to avoid repeated expensive construction of large operators or metadata.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `clebsch_gordan_bypass()`, `clebsch_gordan_general()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- max_rank -maximum L rank for Wigner D functions
- int_ranks -row vector of interaction ranks for which
- product superoperators are required; may
- be empty if only rotation generators are
- needed

## Outputs

- space_basis -lab space basis set descriptor, in
- [L M N] format, giving indices of
- each Wigner function in the basis.
- Lx,Ly,Lz -representations of lab space rotation
- generators in the Wigner function basis,
- to be used in the building of the lab
- space diffusion operator.
- D -a cell array with one element per interac-
- tion rank r in int_ranks, each a cell array
- of Wigner function product superoperators,
- corresponding to multiplication by D[r,M,N]
- of the basis Wigner functions, to be used
- in the building of the spin Hamiltonian
- operator; D{r} has dimensions (2r+1)x(2r+1)
- Automatic caching is implemented -the function would not re-
- compute operator sets that it can find on disk.
- Note: building product superoperators for interaction ranks other
- than 2 calls clebsch_gordan.m, which requires the Java virtual
- machine; cached operator sets load without it.

## Implementation structure

- Wigner D function basis set and rotation generators required by
- the SLE module. Syntax:
- [Lx,Ly,Lz,D,space_basis]=sle_operators(max_rank,int_ranks)
- max_rank -maximum L rank for Wigner D functions
- int_ranks -row vector of interaction ranks for which
- product superoperators are required; may
- be empty if only rotation generators are
- needed
- space_basis - lab space basis set descriptor, in
- [L M N] format, giving indices of
- each Wigner function in the basis.
- Lx,Ly,Lz -representations of lab space rotation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `mfilename()`, `own_path()`, `num2str()`, `exist()`, `load()`, `space_basis()`, `source_states()`, `destin_states()`, `spdiags()`, `clear()`, `destinations()`, `sources()`, `clebsch_gordan_bypass()`, `clebsch_gordan_general()`, `save()`.
