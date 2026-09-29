# kernel/overloads/@rcv/rcv.m

- Signature: `obj=rcv(M)`, `obj=rcv(dim1,dim2)`, or `obj=rcv(R,C,V,dim1,dim2)`
- Source: [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@rcv/rcv.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=rcv/rcv.m)

## Purpose

Constructs an RCV (row-column-value) sparse-matrix object. RCV names the coordinate storage format.

## Representation and dimensions

An RCV object stores parallel column vectors `row`, `col`, and `val`, plus `numRows`, `numCols`, and an `isGPU` flag. The constructor stores row and column indices as `int64`, values as `double`, and dimensions as `int64`.

- `rcv(M)` accepts an existing RCV object unchanged or a numeric matrix. For a matrix it retains both dimensions and obtains the stored entries from `find(M)`.
- `rcv(dim1,dim2)` creates an empty matrix with those dimensions.
- `rcv(R,C,V,dim1,dim2)` stores the supplied coordinate/value entries after columnising the arrays. The three arrays must have equal element counts; `R` and `C` must be finite real integer coordinates inside the matrix dimensions. Dimensions must be finite, real, non-negative integer scalars. `V` must be numeric; the constructor does not impose a real-valued or finite-value check on it.

The explicit-coordinate form retains repeated coordinates as supplied; this constructor neither sorts nor coalesces them. Converting with MATLAB's `sparse` constructor sums contributions at repeated coordinates; see [`sparse`](sparse.md). Coordinates are 1-based matrix row and column indices, not linear indices. The RCV class has no `subsref` element-indexing overload; convert to a MATLAB sparse matrix for ordinary matrix-element indexing.

For products, [`mtimes`](mtimes.md) handles scalar scaling (which retains RCV form) and compatible RCV/MATLAB-sparse matrix products (which return MATLAB sparse form).

## GPU behaviour

A matrix input preserves its CPU/GPU location. In the explicit-coordinate form, if any of `R`, `C`, or `V` is GPU-resident, all three stored entry arrays are placed on the GPU and `isGPU` is true. The dimension-only empty form is CPU-resident.

## Class predicates

The class reports `true` for `isnumeric`, `ismatrix`, and `isfloat`.
