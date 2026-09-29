# kernel/overloads/@ttclass/subsref.m

## Signature

`answer=subsref(ttrain,reference)`

## Behaviour

A dot reference returns one of the stored properties `ncores`, `ntrains`, `sizes`, `ranks`, `coeff`, `cores`, or `tolerance`; an unlisted property or reference type errors. Parentheses indexing accepts exactly two subscripts, for the represented matrix's row and column. Each must be a scalar index; logical scalars are converted to numeric values. Indices must be real positive integers, no greater than `flintmax` or the corresponding matrix dimension. Vector, range, and other advanced indexing are not implemented.

The row and column linear indices are separately expanded into per-core physical coordinates from the last mode back to the first, using the mode sizes. For each tensor-train component, the selected core slices are contracted along their bond indices from the last core towards the first; that scalar contraction is multiplied by the component's coefficient, and the component contributions are added. Thus a parenthesised matrix-element request returns a scalar, not a tensor-train object. A dot request returns the selected property, and any remaining reference levels are applied recursively to that result.

This accessor does not alter ranks, coefficients, or tolerance. It does not round or truncate the representation; it evaluates a single requested element. The method's row/column interpretation follows the tensor train's stored core order and its per-core row and column mode sizes.

## References

- [Source on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/subsref.m)
- [Spin Dynamics Wiki: `ttclass/subsref.m`](https://spindynamics.org/wiki/index.php?title=ttclass/subsref.m)
