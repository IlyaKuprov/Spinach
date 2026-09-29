# kernel/overloads/@cell/minus.m

- Signature: `C=minus(A,B)`

## Behaviour

This overload applies subtraction inside a cell container; it does not subtract or contract the outer cell array as a numeric tensor. With two cells, it requires `isequal(size(A),size(B))` and evaluates `A{n}-B{n}` at each linear index. With one numeric operand, that same operand is applied to each cell entry: `A{n}-B` for cell-minus-numeric, and `A-B{n}` for numeric-minus-cell. The latter order is significant.

There is no cell-level singleton expansion: two cell inputs must have equal sizes. For each pair of contained numeric arrays, dimensional compatibility and any MATLAB implicit expansion come from MATLAB's elementwise minus. This code does not assemble a larger matrix, sum over entries, or impose a sparse-to-full conversion; each contained subtraction is delegated to MATLAB's operator dispatch.

When cells are used to hold orientation-indexed matrices, this overload pairs entries by cell index and subtracts within each entry. It neither represents orientation weights nor combines or contracts across orientations. The function does not special-case `opium` or `polyadic` objects; any inner operation is whatever MATLAB dispatches for those values.

## Inputs and validation

The consistency check accepts each outer operand if it is numeric or a cell array, and checks equal cell-array sizes only when both operands are cells. It does not validate cell contents or the shapes of contained arrays. Invalid inner operations therefore fail through MATLAB's subtraction dispatch. There is no explicit cell broadcasting or shape reconciliation at the outer level.

## Output

A cell array with the same cell shape as the cell operand. The implementation updates a local copy of that operand and returns it as `C`.

## References

- Source: [`kernel/overloads/@cell/minus.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@cell/minus.m)
- Wiki: [cell/minus.m](https://spindynamics.org/wiki/index.php?title=cell/minus.m)
