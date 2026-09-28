# kernel/utilities/reprows.m

- Signature: `B=reprows(A,row_nums,rep_counts)`

## Purpose

Returns a matrix or cell array with selected rows repeated the specified number of times. Rows not named in `row_nums` are kept once. Each entry in `rep_counts` is the total number of copies of its corresponding selected row in the result.

## Physical / mathematical content

This utility performs row selection and replication; it does not apply a physical model or transform numerical values within the rows.

## Numerical / algorithmic content

The function initializes a repetition count of one for every row, replaces counts at the requested indices, expands `1:size(A,1)` with `repelem`, and indexes `A` with that expanded row list. Selected row indices must be unique, positive integers within the row dimension; repetition counts must be matching positive integers.

## Parameters / inputs

- `A` - numeric matrix or cell array.
- `row_nums` - vector of distinct positive integer row indices, no larger than `size(A,1)`.
- `rep_counts` - vector with one positive integer per selected row, specifying its total number of copies.

## Outputs

- `B` - same type as `A`, with the rows expanded according to `row_nums` and `rep_counts`.

## Implementation structure

Input checks reject other types for `A`, non-real or non-finite or non-integer row indices, duplicate or out-of-range row indices, and repetition-count vectors that do not match `row_nums` or contain invalid counts. The expanded index vector is then used for the row extraction.

## Reference

[Spin Dynamics Wiki: reprows.m](https://spindynamics.org/wiki/index.php?title=reprows.m)
