# kernel/overloads/@cell/mtimes.m

- Signature: `C=mtimes(A,B)`

## Behaviour

Exactly one operand must be a cell array and the other numeric. For cell-left multiplication, the operation is `A{n}*B` for every entry; for numeric-left multiplication, it is `A*B{n}`. Operand order is preserved, so these cases are not interchangeable. Both-cell inputs are rejected by the overload.

The product at each entry is MATLAB matrix multiplication, not elementwise multiplication. Its inner matrix dimensions must satisfy MATLAB's multiplication rules; this wrapper does not broadcast or contract across the cell dimensions. If cells hold orientation-indexed matrices, the same numeric matrix is multiplied into each entry independently, and no reduction or orientation weighting is performed. Sparse operands are passed to the contained `*` operations without an explicit conversion; result storage follows MATLAB's operation. The implementation does not define separate `opium` or `polyadic` semantics: if an inner value invokes another class overload, that behaviour is outside this wrapper.

## Inputs and validation

The consistency check accepts outer operands that are numeric or cells. The branch logic permits a cell paired with a numeric operand and errors for both cells; cell contents are not separately validated. MATLAB's multiplication dispatch determines whether each contained product is defined and dimensionally valid. No cell-array broadcasting rule is implemented.

## Output

A cell array with the shape of the cell operand, with each entry replaced by its ordered matrix product.

## References

- Source: [`kernel/overloads/@cell/mtimes.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@cell/mtimes.m)
- Wiki: [cell/mtimes.m](https://spindynamics.org/wiki/index.php?title=cell/mtimes.m)
