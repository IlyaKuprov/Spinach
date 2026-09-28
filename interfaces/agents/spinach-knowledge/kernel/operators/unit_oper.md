# kernel/operators/unit_oper.m

- Signature: `A=unit_oper(spin_system)`

## Purpose

Returns a sparse identity operator with dimensions appropriate to the current formalism and basis.

## Physical / mathematical content

- In `sphten-liouv`, the dimension is the number of rows in `spin_system.bas.basis`.
- In `zeeman-hilb` and `zeeman-wavef`, the dimension is the product of the spin multiplicities.
- In `zeeman-liouv`, the dimension is the square of that product.

## Numerical / algorithmic content

- Constructs the sparse identity matrix with `speye`.

## Parameters / inputs

- `spin_system` — Spinach data object containing basis information; call `basis.m` first.

## Outputs

- `A` — Sparse identity matrix of the appropriate dimension.

## Implementation structure

- Checks that `spin_system.bas.formalism` is present, then selects the dimension by formalism. An unknown formalism raises an error.

## Reference

- https://spindynamics.org/wiki/index.php?title=unit_oper.m