# kernel/overloads/@cell/blkdiag.m

- Signature: `C=blkdiag(A,B)`

This overload takes two cell matrices and returns a cell matrix whose row and column counts are the sums of the corresponding input dimensions. The cells of `A` occupy the upper-left block and those of `B` the lower-right block; the off-diagonal positions remain empty cells. Cell contents are placed, not added or combined, so this is block placement on the cell grid rather than numeric block-diagonal assembly of the contents.

The checks require both inputs to be cell arrays and matrices. The implementation does not check or transform the contents and defines no broadcasting rule; the block sizes follow the input cell-array dimensions.

## References

- MATLAB source: [`kernel/overloads/@cell/blkdiag.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@cell/blkdiag.m)
- Spinach Wiki: [`cell/blkdiag.m`](https://spindynamics.org/wiki/index.php?title=cell/blkdiag.m)
