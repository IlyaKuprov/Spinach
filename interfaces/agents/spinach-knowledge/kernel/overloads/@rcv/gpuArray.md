# kernel/overloads/@rcv/gpuArray.m

- Signature: `obj=gpuArray(obj)`

## Purpose

Moves an RCV sparse matrix's stored arrays to GPU memory when it is not already GPU-resident.

## Physical / mathematical content

The transfer changes where the row indices, column indices, and values are stored, not the represented sparse matrix.

## Numerical / algorithmic content

If obj.isGPU is false, gpuArray is applied to obj.row, obj.col, and obj.val, and the flag is set to true. A GPU-resident input is left unchanged.

## Parameters / inputs

- obj -an RCV sparse matrix

## Outputs

- obj -the same matrix with data stored on GPU

## Implementation structure

- Requires obj to be an RCV object.
- Transfers the row, column, and value arrays only when obj.isGPU is false.
- Sets obj.isGPU after transferring those arrays.
