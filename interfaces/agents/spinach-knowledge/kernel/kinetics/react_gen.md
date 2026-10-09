# kernel/kinetics/react_gen.m

## Product-row reaction compilation

`maps=react_gen(spin_system,reaction)` accepts a validated reaction record and the compiled spherical-tensor direct sum. It returns one index matrix per product occurrence: the first column is a global destination row, and the remaining columns are global source rows, one per reactant occurrence.

Compilation enumerates product rows, pulls their descriptors back through the atom matching, and looks each source descriptor up in its local substance basis. Unmatched source spins are traced out; unmatched product spins can only arrive at identity. The unit row is retained, making concentrations part of the same reaction map as spin orders. A spin-free reactant contributes its sole unit coordinate. Repeated product indices produce repeated maps with the required stoichiometry.

The compiler never forms a Cartesian product of reactant bases. Product rows excluded because a source descriptor was truncated are counted and reported separately from orders on unmatched product spins. `kinetics` applies rates, drains, selectors, and the additive or product closure to these index lists.

The two-column global matching cannot identify different molecular occurrences of a repeated reactant. A repeated reactant whose spins appear in matching therefore raises `Spinach:react_gen:repeatedMatching`; repeated spin-free or wholly traced reactants remain supported. Matched repeated spin-bearing products likewise raise `Spinach:react_gen:repeatedProductMatching`, because global destination labels do not identify molecular occurrences; repeated unlabelled products retain their stoichiometric multiplicity. Occurrence-resolved matching is not guessed.
