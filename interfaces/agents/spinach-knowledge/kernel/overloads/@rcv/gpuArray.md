# kernel/overloads/@rcv/gpuArray.m

[GitHub source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@rcv/gpuArray.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=rcv/gpuArray.m)

- Signature: `obj=gpuArray(obj)`

## Purpose

Moves a CPU-resident RCV sparse matrix's stored arrays to GPU memory.

## Storage and behaviour

RCV stores row indices, column indices, and corresponding values in parallel arrays; `numRows` and `numCols` retain the matrix shape. The overload checks that `obj` is an `rcv` object. If `obj.isGPU` is false, it eagerly applies `gpuArray` to `obj.row`, `obj.col`, and `obj.val`, then sets `obj.isGPU` to true. If the flag is already true, it leaves the object unchanged. The row and column counts are not reassigned, so the represented dimensions remain unchanged. This transfers the stored arrays; it does not build a MATLAB sparse or dense matrix.

The values are transferred without conjugation or scalar expansion/broadcasting.

## Input

- `obj` - an RCV sparse matrix. The explicit check is object type only.

## Output

- `obj` - the same RCV matrix with its stored arrays on the GPU when it was CPU-resident.
