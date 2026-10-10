# kernel/summaries/summary_basis.m

Source: [kernel/summaries/summary_basis.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/summaries/summary_basis.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=summary_basis.m)

- Signature: `summary_basis(spin_system)`

## Purpose

Reports the final basis-set summary for a Spinach system. For each listed basis state, each spin's stored basis entry is converted by `lin2lm` to its irreducible spherical-tensor quantum-number pair `(L,M)`; the table columns are the spin indices.

## Numerical content

Each substance is reported separately using `bas.nstates(n)` and its local descriptor `bas.basis{n}`. Column headings identify global spins in `chem.parts{n}`, while state numbers include the compiled offset. The percentage uses the full local Hilbert dimension for Hilbert/wavefunction formalisms and its square for Liouville formalisms. Non-sphten formalisms print dimensions without tensor labels. If `nstates` is greater than `spin_system.tols.basis_hush`, the detailed state table is suppressed; the final dimension and percentage are still reported. The labels and counts have no physical units.

## Parameters / inputs

- `spin_system` - Spinach spin system structure.

## Output and side effects

Writes the heading, optional per-state `(L,M)` table, and final dimension/percentage through `report.m`, which routes the text to the console or the configured output. The function checks that `spin_system` is a structure and otherwise does not modify it.

Legacy global `bas.basis` matrices and `bas.irrep` fields are rejected at this entry point with named errors pointing to per-substance `bas.basis{n}`/`bas.offsets` and `bas.sym_fact(n)` symmetry data. Compiled structures remain ordinary MATLAB structs; arbitrary external dot reads are not intercepted.
