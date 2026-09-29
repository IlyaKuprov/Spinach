# kernel/overloads/@cell/plus.m

- Signature: `C=plus(A,B)`

## Behaviour

This overload adds values within a cell container; it does not add or reduce the outer cell array as a numeric tensor. With two cells, it requires `isequal(size(A),size(B))` and evaluates `A{n}+B{n}` at each linear index. With one numeric operand, that same operand is added to each cell entry: `A{n}+B` for cell-plus-numeric, and `A+B{n}` for numeric-plus-cell.

There is no cell-level singleton expansion: two cell inputs must have equal sizes. For each pair of contained arrays, dimensional compatibility and any MATLAB implicit expansion come from MATLAB's elementwise plus. This implementation does not assemble a larger matrix or sum across cells, and makes no sparse/full conversion; storage behaviour comes from each underlying MATLAB operation.

For orientation-indexed cell data, addition pairs entries by cell index and acts independently within each entry. It does not apply orientation weights or reduce across orientations. The wrapper does not define special `opium` or `polyadic` operations; inner operations follow MATLAB's dispatch for the contained values.

## Inputs and validation

Each outer operand must be numeric or a cell array. Equal cell sizes are checked only when both operands are cells. The contents and their array dimensions are not checked by the wrapper; MATLAB's addition dispatch determines whether each inner operation is valid. No broadcasting between differently shaped cell arrays is implemented.

## Output

A cell array with the shape of the cell operand. The implementation updates a local copy of that operand and returns it as `C`.

## References

- Source: [`kernel/overloads/@cell/plus.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@cell/plus.m)
- Wiki: [cell/plus.m](https://spindynamics.org/wiki/index.php?title=cell/plus.m)
