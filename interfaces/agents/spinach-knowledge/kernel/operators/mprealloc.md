# kernel/operators/mprealloc.m

- Signature: `A=mprealloc(spin_system,nnzpc)`
- Direct source: [kernel/operators/mprealloc.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/operators/mprealloc.m)
- Wiki: [mprealloc.m](https://spindynamics.org/wiki/index.php?title=mprealloc.m)

## Purpose

Allocates an all-zero sparse square matrix sized for the active Spinach formalism, reserving an estimated number of nonzeros per column. This routine allocates storage only: it does not construct matrix elements, define an operator's action, or propagate a state.

## Dimension and formalism

The function reads `spin_system.bas.formalism` and chooses the square dimension as follows:

- `sphten-liouv`: `size(spin_system.bas.basis,1)`.
- `zeeman-wavef` and `zeeman-hilb`: `prod(spin_system.comp.mults)`.
- `zeeman-liouv`: `prod(spin_system.comp.mults.^2)`.

For the first case the dimension is the row count of the supplied basis matrix; in the other cases it is calculated from the spin multiplicities. The routine does not create or reorder any basis, so the basis ordering is whatever the selected formalism already uses.

## Allocation and input checks

For the selected dimension `d`, the implementation calls `spalloc(d,d,nnzpc*d)`; the third argument is reserved sparse storage, not a count of nonzeros already present. The function checks that `spin_system.bas.formalism` exists and that `nnzpc` is numeric, real, scalar, and integer-valued before allocation. An unrecognised formalism raises an error. The documentation describes `nnzpc` as the expected nonzero count per column.
