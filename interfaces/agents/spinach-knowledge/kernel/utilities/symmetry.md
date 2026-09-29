# kernel/utilities/symmetry.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/symmetry.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/symmetry.m)

## Purpose

Permutation symmetry treatment. Compiles character tables of composite symmetry groups, builds the permutation table for each spin state in the basis, and builds projectors into the irreducible representations of the direct product symmetry group.

## Behaviour

- Syntax: `spin_system=symmetry(spin_system,bas)`.
- This is a service function of the Spinach kernel that should not be called directly; it is called by `basis.m`.
- Non-Abelian groups and multi-dimensional irreps are supported; edit `perm_group.m` to add your own groups.
- Consistency is enforced by an internal `grumble` subfunction which validates the symmetry parameters in `bas` (see Inputs and outputs).
- If `'symmetry'` is listed in `spin_system.sys.disable`, the function issues a warning that symmetry factorisation is disabled and writes empty cells to `spin_system.comp.sym_group`, `spin_system.comp.sym_spins`, and sets `spin_system.comp.sym_a1g_only` to true.
- Otherwise, the fields `sym_group`, `sym_spins` and `sym_a1g_only` are taken from `bas` when present. When `sym_a1g_only` is not given in `bas`, it defaults to false for the `zeeman-hilb` and `zeeman-wavef` formalisms and true otherwise.
- If symmetry groups are declared, a permutation symmetry summary is printed via `summary_symmetry`, and `validate_sym` checks that the interactions respect the declared symmetry.
- If more than one symmetry group is declared, the constituent groups are lifted from the `perm_group` database, the direct product character table is computed as a Kronecker product of the individual character tables, and the direct product element list and group order are computed. The spin lists from `sym_spins` are concatenated.
- If exactly one symmetry group is declared, it is lifted from the `perm_group` database and the corresponding spin list from `sym_spins` is used.
- If no symmetry group is declared, the function reports that no symmetry information is available.
- When a group is available, the permutation table over all group operations is computed with a `parfor` loop; for the `zeeman-liouv` formalism the permutation is extended to include the adjoint spin indices (`group_element+nspins`). Row indexing of the permuted basis uses `spsortrows`.
- In the fully symmetric (`sym_a1g_only` true) mode, only the A1g irrep is retained: the permutation table is pruned to unique sorted rows, a coefficient matrix is built, normalised, and returned as the projector with its dimension.
- In the full symmetry treatment mode, the function loops over all irreducible representations. For each irrep it builds a sparse transformation matrix from the permutation table and the character table (skipping zero characters), removes sign ambiguity, cleans up using `spin_system.tols.liouv_zero`, removes zero and identical columns (`spunicols`), and then either orthogonalises multi-dimensional irreps (using `scomponents` to find non-orthogonal subspaces and `orth` to orthogonalise them) or normalises the SALCs for one-dimensional irreps. Zero-dimensional irreps are removed at the end.
- Progress and diagnostics are reported to the user via `report`, including the number of irreps, their dimensions, the number of symmetry operations in the group direct product, and the number of states per irrep.

## Inputs and outputs

**Inputs**

- `spin_system` — Spinach spin system description object produced as described in the spin system and basis specification sections of the online manual.
- `bas` — basis input structure described in the basis specification section of the manual. Recognised symmetry-related fields:
  - `bas.sym_group` — cell array of strings naming the symmetry groups; supported group names are `'S2'`, `'S3'`, `'S4'`, `'S4A'`, `'S5'`, `'S6'`, `'S6A'`, `'S8A'`.
  - `bas.sym_spins` — cell array of numeric vectors of spin indices, one per symmetry group; must have the same number of elements as `bas.sym_group`. Each vector must contain at least two integer spin labels within `1..nspins`, the vectors must not intersect, and each group must not cross chemical substance boundaries (`spin_system.chem.parts`).
  - `bas.sym_a1g_only` — irrep composition switch; allowed values are 0 and 1. It must be specified alongside `bas.sym_group`.
  - `bas.sym_spins` and `bas.sym_a1g_only` may not be specified without `bas.sym_group`.

**Outputs**

- `spin_system.bas.irrep(n).projector` — projector matrices into each irreducible representation.
- `spin_system.bas.irrep(n).dimension` — dimension of each irreducible representation.
- The function also populates `spin_system.comp.sym_group`, `spin_system.comp.sym_spins` and `spin_system.comp.sym_a1g_only`.

## References

- Spinach online manual: [https://spindynamics.org/wiki/index.php?title=symmetry.m](https://spindynamics.org/wiki/index.php?title=symmetry.m)
