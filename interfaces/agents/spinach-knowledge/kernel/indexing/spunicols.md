# kernel/indexing/spunicols.m

Source: [MATLAB implementation](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/indexing/spunicols.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=spunicols.m).

- Signature: `A=spunicols(A)`

## Purpose

Returns one copy of every distinct column of a sparse matrix. This is a column-set operation, not a physical model and not a line-shape or chemical-flow calculation.

## Mathematical mapping

For an input `A` with `m` rows and `n` columns, regard each column `a_j` as an element of `R^m`. Define `a_i ~ a_j` exactly when `a_i=a_j`. The result `B` contains one column for each equivalence class, so its shape is `m × u`, where `u` is the number of distinct input columns. In code the mapping is `B = unique(A.','rows').'`.

MATLAB's `unique(...,'rows')` uses its default sorted ordering; the output columns are therefore ordered by the corresponding sorted rows, not kept in first-occurrence order. The function has one output only: it does not return the source-column indices or an index map. It applies no tolerance, normalisation, or physical-unit conversion.

## Input and output

The source contract is a sparse, real, double, two-dimensional matrix. The MATLAB fallback returns a sparse real double matrix made from the unique columns. The entries' interpretation and units are inherited unchanged from the caller; the operation itself is unit-agnostic.

## Source guard

The local `grumble` check rejects inputs that are not numeric, sparse, real, class `double`, and a matrix, with the error “A must be a sparse real double matrix.” The source is documented as the MATLAB fallback for the compiled MEX implementation.
