# kernel/overloads/@rcv/gather.m

[GitHub source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@rcv/gather.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=rcv/gather.m)

- Signature: `A=gather(A)`

## Purpose

Moves a GPU-resident RCV sparse matrix's stored arrays to CPU memory.

## Storage and behaviour

RCV stores row indices, column indices, and corresponding values in parallel arrays; `numRows` and `numCols` retain the matrix shape. The overload checks that `A` is an `rcv` object. If `A.isGPU` is true, it eagerly applies `gather` to `A.row`, `A.col`, and `A.val`, then sets `A.isGPU` to false. If the flag is already false, it leaves the object unchanged. The row and column counts are not reassigned, so the represented dimensions remain unchanged. This transfers the stored arrays; it does not build a MATLAB sparse or dense matrix.

The values are transferred without conjugation or scalar expansion/broadcasting.

## Input

- `A` - an RCV sparse matrix. The explicit check is object type only.

## Output

- `A` - the same RCV matrix with its stored arrays on the CPU when it was GPU-resident.
