# kernel/overloads/@cell/times.m

- MATLAB implementation: [kernel/overloads/@cell/times.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@cell/times.m)

- Signature: `C=times(A,B)`
- Reference: [Spinach Wiki: `cell/times.m`](https://spindynamics.org/wiki/index.php?title=cell/times.m)

## Purpose and supported operands

This overload applies MATLAB multiplication across the entries of a cell array. Exactly one argument must be a cell array and the other a numeric array; two cell arrays are rejected. A non-cell argument must be numeric. The result `C` is a cell array with the same shape as the cell operand.

## Pairing and multiplication rules

Use a numeric scalar to apply the same factor to every cell. Alternatively, the numeric array must have exactly the same dimensions as the cell array—not merely the same number of elements—and its entries are paired by linear index with the corresponding cells. This works in either operand order:

- Cell on the left: each entry becomes `A{n}*B` for scalar `B`, or `A{n}*B(n)` for a matching numeric array.
- Cell on the right: each entry becomes `A*B{n}` for scalar `A`, or `A(n)*B{n}` for a matching numeric array.

The operator inside each cell is MATLAB matrix multiplication (`*`), not element-wise multiplication of the cell contents. Consequently each cell value must support the corresponding matrix product. The implementation checks operand types, rejects cell/cell input, and errors when a non-scalar numeric array's size differs from the cell array's size; it does not broadcast non-scalar arrays or transpose/reorder them to make shapes match.
