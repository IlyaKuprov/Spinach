# kernel/kinetics/react_gen.m

- Signature: `G=react_gen(spin_system,reaction)`
- Direct MATLAB source: [`kernel/kinetics/react_gen.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/kinetics/react_gen.m)
- Existing Wiki: [`react_gen.m`](https://spindynamics.org/wiki/index.php?title=react_gen.m)

## Purpose and inputs

Builds state-space reaction generators for the listed reactant and product parts of `spin_system.chem.parts`. `reaction.reactants` and `reaction.products` are integer part indices. Each row of the two-column `reaction.matching` pairs a reactant spin index (column 1) with its corresponding product spin index (column 2).

The guard rejects overlapping reactant/product part lists, unequal total numbers of spins on the two sides, or a matching map whose left/right indices fail to cover the respective reactant/product spin lists. During construction, each populated basis state must resolve to exactly one host substance; a destination state is looked up in the existing basis.

## State-index map

For each basis row `n`, the code identifies the host part from the spins active in that state. If that part is a reactant, it records a drain `-1` at (n,n) in that reactant's generator. It forms the product state by copying entries of the source basis row from `reaction.matching(:,1)` to `reaction.matching(:,2)`, then searches for that complete row in `spin_system.bas.basis`. If found at index `d`, it records a fill `+1` at (d,n); if absent, no fill entry is added. Thus for a mapped state `n→d`, the column-indexed action is `G_j[:,n] = -e_n + e_d`; when the destination row is not present, it is `G_j[:,n] = -e_n`.

## Output and lifecycle

`G` is a cell array with one matrix per listed reactant. Each member is an `nstates×nstates` complex sparse matrix, where `nstates=size(spin_system.bas.basis,1)`. The stored coefficients are dimensionless `-1/+1` generator entries; a reaction rate is not an input to this function and must be applied by the caller. There is no normalisation step. The routine collects drain/fill indices in parallel, removes unused rows, and assembles each reactant's sparse matrix.

## Source Wiki

[Spinach Wiki: `react_gen.m`](https://spindynamics.org/wiki/index.php?title=react_gen.m)
