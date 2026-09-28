# kernel/overloads/@rcv/gather.m

- Signature: `A=gather(A)`

## Purpose

Moves an RCV sparse matrix's stored arrays from GPU to CPU memory when it is GPU-resident.

## Physical / mathematical content

Gathering changes the location of the row indices, column indices, and values, not the represented sparse matrix.

## Numerical / algorithmic content

If A.isGPU is true, gather is applied to A.row, A.col, and A.val, and the flag is set to false. A CPU-resident input is left unchanged.

## Parameters / inputs

- A -an RCV sparse matrix

## Outputs

- A -the same matrix with data stored on the CPU

## Implementation structure

- Requires A to be an RCV object.
- Gathers the row, column, and value arrays only when A.isGPU is true.
- Clears A.isGPU after transferring those arrays.
