# kernel/utilities/kronm.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/kronm.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=kronm.m)

## Purpose and interface

`x=kronm(Q,x)` applies `kron(Q{1},kron(Q{2},...))*x` without forming the Kronecker product. `Q` is a cell array of matrix factors; `x` is a numeric vector or matrix whose row count equals the product of the factor column counts. The output has the product of their row counts and the same number of columns as the input. Rectangular factors are supported.

The contraction reshapes each tensor dimension into matrix rows, applies the corresponding factor to all remaining columns, and restores the tensor layout. Factors are applied in reverse order to respect MATLAB column-major Kronecker ordering. `opium` factors use their scalar coefficient without forming identity matrices.

## Implicit factors

A factor may instead be a construction description with `action` and `dims` fields, as supplied to `polyadic`. `dims=[nrows ncols]` provides the tensor dimensions, and `action(block)` applies the factor to a numeric matrix with `ncols` rows. The existing polyadic `mtimes` overload supplies the action and dimensions internally; an FFT or other action needs no string-command protocol. The contraction itself does not need the adjoint field.

The function validates the cell structure, scalar implicit descriptions with function-handle actions, the factor dimensions, and the numeric right-hand side. The caller is responsible for matching dimensions and supplying linear actions with the documented output shape. No periodic-grid or spin-system assumptions are made here.
