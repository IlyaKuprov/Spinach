# kernel/summaries/summary_basis.m

Source: [kernel/summaries/summary_basis.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/summaries/summary_basis.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=summary_basis.m)

- Signature: `summary_basis(spin_system)`

## Purpose

Reports the final basis-set summary for a Spinach system. For each listed basis state, each spin's stored basis entry is converted by `lin2lm` to its irreducible spherical-tensor quantum-number pair `(L,M)`; the table columns are the spin indices.

## Numerical content

The reported dimension is the number of rows in `spin_system.bas.basis`. The final percentage is `100 * nstates / (prod(spin_system.comp.mults)^2)`, printed as a percentage of the full state space. If `nstates` is greater than `spin_system.tols.basis_hush`, the detailed state table is suppressed; the final dimension and percentage are still reported. The labels and counts have no physical units.

## Parameters / inputs

- `spin_system` - Spinach spin system structure.

## Output and side effects

Writes the heading, optional per-state `(L,M)` table, and final dimension/percentage through `report.m`, which routes the text to the console or the configured output. The function checks that `spin_system` is a structure and otherwise does not modify it.
