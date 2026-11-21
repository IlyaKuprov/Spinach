# kernel/operators/superop.m

- Source: [kernel/operators/superop.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/operators/superop.m)
- Wiki: [superop.m on the Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=superop.m)
- Signature: `A=superop(spin_system,opspec,side)`

## Purpose and basis

Builds the matrix representation of left or right multiplication of a density operator in the selected `sphten-liouv` basis. It requires basis information in `spin_system.bas`; the represented matrix acts on the substance-local descriptors in `spin_system.bas.basis{n}`, rather than directly on the Hilbert-space state vector. Active global spins must belong to one substance; `which_subst` rejects a cross-substance product. Columns are translated through `chem.parts{n}`, and the triplet row and column indices are shifted by `bas.offsets(n)` into the direct sum.

## Side and action

- `side='left'` selects the left-product action; `side='right'` selects the right-product action.
- `side='comm'` combines the left action with the negative of the right action, i.e. the commutator action. Rank-exchange parity removes cancelling product paths before a single basis-matching traversal.
- `side='acomm'` retains even-total-parity product paths and doubles their coefficients, giving the anticommutator action.

The two commutator-specific branches discard local transitions when the sum of the active-spin source indices or the sum of their destination indices is zero. The general left/right branches do not apply that filter. With an all-zero `opspec`, left/right return the identity, the commutator returns zero, and the anticommutator returns twice the identity.

## Specification, coefficients, and output

`opspec` is a row vector with one operator-state index per spin. For each active spin, the implementation selects the multiplicity-specific left-product table `spin_system.bas.lpst` or right-product table `spin_system.bas.rpst`; each table stores `[n m k value]` coordinate records. The `lp_offsets` and `rp_offsets` arrays locate contiguous records for each first tensor index, avoiding full-table scans. Selecting `n=opspec+1` supplies the two remaining indices and coefficients directly, without expanding a dense page or calling `find` on it. Tensor indices are then converted to Spinach's zero-based state labels. For multiple active spins, the indices and coefficients are combined by Kronecker products, then matched against rows of the hosting `spin_system.bas.basis{n}` to assemble the action in the selected basis. No extra normalisation factor is applied in this assembly.

- `opspec`: Spinach operator specification described in Sections 2.1 and 3.3 of [the cited paper](http://dx.doi.org/10.1016/j.jmr.2010.11.008). It must contain one integer for each spin, with operator-state indices permitted by the corresponding multiplicity.
- `side`: `'left'`, `'right'`, `'comm'`, or `'acomm'`.
- `A`: three-column XYZ sparse triplets: row index, column index, and value. This is not MATLAB's CSC storage format; if the selected transitions are empty, the source returns the zero triplet `[1 1 0]`.

Direct calls to this general function are not usually required; the source recommends the friendlier `operator()` function.

For unequal source and destination subspaces, active labels in the full source rows are replaced by the target transition labels, then `hash_ismember` matches these rows against the selected cached `bas.row_hash{n}` destination rows. This retains the spectator-row equality condition without hashing references again. Equal subspaces still use the direct-index shortcut.

For product commutators, odd parity is required for the sum of operator, source, and destination ranks over all active spins; anticommutators retain even parity. This is not a product of single-spin commutators. The retained left-product coefficients are doubled before basis matching.
