# kernel/operators/sin_tran.m

- Signature: `A=sin_tran(dim)`

## Purpose

Returns single-transition operators spanning the space of matrices of dimension `dim`. Each operator is a sparse matrix with one nonzero entry; the cell array uses serpentine indexing.

## Physical / mathematical content

For a 4-by-4 matrix, cell-array positions map to matrix entries as follows:

```text
 1   3   6  10
 2   5   9  13
 4   8  12  15
 7  11  14  16
```

The pattern continues for larger matrices.

## Numerical / algorithmic content

The function allocates a `dim^2`-element cell array and constructs each operator as a sparse matrix with one unit entry. Operators are converted to complex type to avoid expensive reallocations later. Construction uses `parfor`.

## Parameters / inputs

- `dim` — matrix dimension; must be a positive real integer.

## Outputs

- `A` — cell array of `dim^2` complex sparse matrices, ordered as shown above.

## Reference

<https://spindynamics.org/wiki/index.php?title=sin_tran.m>