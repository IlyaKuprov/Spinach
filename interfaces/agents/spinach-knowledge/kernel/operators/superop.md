# kernel/operators/superop.m

- Source: [kernel/operators/superop.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/operators/superop.m)
- Wiki: [superop.m on the Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=superop.m)
- Signature: `A=superop(spin_system,opspec,side)`

## Purpose and basis

Builds the matrix representation of left or right multiplication of a density operator in the selected `sphten-liouv` basis. It requires basis information in `spin_system.bas`; the represented matrix acts on the basis stored in `spin_system.bas.basis`, rather than directly on the Hilbert-space state vector. If that basis has `N` rows, the action is an `N`-by-`N` superoperator returned as sparse triplets, not as a dense matrix.

## Side and action

- `side='left'` selects the left-product action; `side='right'` selects the right-product action.
- `side='comm'` combines the left action with the negative of the right action, i.e. the commutator action. The implementation obtains these from the `leftofcomm` and `rightofcomm` branches and negates the latter's values.
- `side='acomm'` concatenates the left and right actions, i.e. the anticommutator action.

The two commutator-specific branches discard local transitions when the sum of the active-spin source indices or the sum of their destination indices is zero. The general left/right branches do not apply that filter. With an all-zero `opspec`, the function takes the identity-operator shortcut.

## Specification, coefficients, and output

`opspec` is a row vector with one operator-state index per spin. For each active spin, the implementation selects the multiplicity-specific left-product table `spin_system.bas.lpst` or right-product table `spin_system.bas.rpst`; table row/column indices are converted from Spinach's zero-based operator-state indexing, and table values supply the local structure coefficients. For multiple active spins, the indices and coefficients are combined by Kronecker products, then matched against rows of `spin_system.bas.basis` to assemble the action in the selected basis. No extra normalisation factor is applied in this assembly.

- `opspec`: Spinach operator specification described in Sections 2.1 and 3.3 of [the cited paper](http://dx.doi.org/10.1016/j.jmr.2010.11.008). It must contain one integer for each spin, with operator-state indices permitted by the corresponding multiplicity.
- `side`: `'left'`, `'right'`, `'comm'`, or `'acomm'`.
- `A`: three-column XYZ sparse triplets: row index, column index, and value. This is not MATLAB's CSC storage format; if the selected transitions are empty, the source returns the zero triplet `[1 1 0]`.

Direct calls to this general function are not usually required; the source recommends the friendlier `operator()` function.
