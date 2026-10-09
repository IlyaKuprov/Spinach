# kernel/coherence.m

- Signature: `rho = coherence(spin_system,rho,spec)`
- Implementation: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/coherence.m>

## Contract and supported bases

Keeps only the specified coherence orders in a state vector, a horizontal stack of state vectors, or (for `zeeman-hilb`) a density matrix or horizontal stack of density matrices. Supported formalisms are `sphten-liouv`, `zeeman-liouv`, and `zeeman-hilb`; the source also notes support for Fokker-Planck direct products in the Liouville-space formalisms.

Let `N = spin_system.bas.offsets(end)`. The implementation folds the input into a two-dimensional working array with `N` rows in the Liouville formalisms and `N^2` rows in `zeeman-hilb`; the remaining dimension is `numel(rho)` divided by that row count. It applies a row mask to all columns, then reshapes the result back to the input's original dimensions. For `zeeman-hilb`, the row count `N^2` represents the stretched density-matrix elements.

Coherence orders are obtained from the basis projections: `lin2lm` supplies projection orders from each local `bas.basis{n}` for `sphten-liouv`, with selected global spins mapped through `chem.parts{n}`. Both Zeeman cases construct exact half-integer projection labels from the spin multiplicities and form ket-minus-bra differences. Multi-substance Zeeman selection is explicitly rejected pending its formalism implementation. For each specification, the function sums those orders over the selected spins and retains rows whose sum belongs to the supplied order vector. Masks from separate specifications are intersected, so every specification must match.

## Specification

`spec` is a cell array of nested pairs: a spin selector and a vector of real integer coherence orders. A selector can be an isotope string in the system, `'electrons'`, `'nuclei'`, `'all'`, or a vector of spin numbers. For example, `{{'13C',[1 -1]},{'1H',-1}}` keeps entries with order 1 or -1 on 13C and order -1 on 1H. The selected spin sets and order lists are evaluated independently for each pair, then ANDed by mask intersection.

The state array must be numeric. The implementation warns through `report` if the filtered result has 1-norm below `1e-10`, with the message that all magnetisation appears to have been destroyed.

## Reference

- [Spinach Wiki: `coherence.m`](https://spindynamics.org/wiki/index.php?title=coherence.m)

Legacy global `bas.basis` matrices and `bas.irrep` fields are rejected at this entry point with named errors pointing to per-substance `bas.basis{n}`/`bas.offsets` and `bas.sym_fact(n)` symmetry data. Compiled structures remain ordinary MATLAB structs; arbitrary external dot reads are not intercepted.
