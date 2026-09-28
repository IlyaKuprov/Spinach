# kernel/utilities/symmetry.m

- Signature: `spin_system=symmetry(spin_system,bas)`

## Purpose

Treats permutation symmetry in a Spinach basis. It obtains character tables for the specified symmetry groups, constructs their direct product when necessary, builds a table of how group operations permute basis states, and constructs projectors into irreducible representations.

## Parameters / inputs

- `spin_system`: Spinach spin-system description object, prepared as described in the spin-system and basis-specification sections of the online manual.
- `bas`: Basis input structure described in the basis-specification section of the manual. Optional symmetry settings include `bas.sym_group` (a cell array of group-name strings), `bas.sym_spins` (a corresponding cell array of spin-index vectors), and `bas.sym_a1g_only` (whether to retain only the fully symmetric representation). Supported group names are `S2`, `S3`, `S4`, `S4A`, `S5`, `S6`, `S6A`, and `S8A`.

## Outputs

- `spin_system.bas.irrep(n).projector`: Projector matrix for an irreducible representation.
- `spin_system.bas.irrep(n).dimension`: Number of states in that representation's projected basis.

## Numerical / algorithmic content

- When multiple symmetry groups are specified, their character tables and permutation elements are combined as a direct product.
- Group operations permute the spin labels of basis states; in `zeeman-liouv` formalism, the permutation is applied to both halves of the basis state.
- With `sym_a1g_only` enabled, the function groups symmetry-related basis states and constructs a normalized projector for the fully symmetric representation. Otherwise it constructs symmetry-adapted linear combinations for all irreducible representations, orthogonalizes overlapping subspaces of multidimensional representations, and removes zero-dimensional results.
- Unless explicitly supplied, `sym_a1g_only` defaults to false for `zeeman-hilb` and `zeeman-wavef`, and true for other formalisms. The symmetry treatment can be disabled through `spin_system.sys.disable`.

This is a Spinach kernel service function called by `basis.m`, not intended for direct use. Non-Abelian groups and multidimensional irreducible representations are supported; edit `perm_group.m` to add groups.

https://spindynamics.org/wiki/index.php?title=symmetry.m