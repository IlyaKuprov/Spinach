# kernel/operators/mprealloc.m

- Signature: `A=mprealloc(spin_system,nnzpc)`
- Direct source: [kernel/operators/mprealloc.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/operators/mprealloc.m)
- Wiki: [mprealloc.m](https://spindynamics.org/wiki/index.php?title=mprealloc.m)

## Purpose

Allocates an all-zero sparse square matrix sized for the active Spinach formalism, reserving an estimated number of nonzeros per column. This routine allocates storage only: it does not construct matrix elements, define an operator's action, or propagate a state.

## Dimension and formalism

The dimension is `spin_system.bas.offsets(end)`, compiled by `basis` as the sum of the substance dimensions in the selected formalism. No basis is constructed or reordered.

## Allocation and input checks

For the selected dimension `d`, the implementation calls `spalloc(d,d,nnzpc*d)`; the third argument is reserved sparse storage, not a count of nonzeros already present. The function checks that `spin_system.bas.formalism` exists and that `nnzpc` is numeric, real, scalar, and integer-valued before allocation. The documentation describes `nnzpc` as the expected nonzero count per column.
