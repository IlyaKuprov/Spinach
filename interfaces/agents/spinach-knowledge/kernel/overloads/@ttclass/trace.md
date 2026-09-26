# kernel/overloads/@ttclass/trace.m

- Signature: `tttrace=trace(tt)`

## Purpose

Computes the trace of a tensor train operator. Syntax: tttrace=trace(tt)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

## Parameters / inputs

- tt -tensor train operator

## Outputs

- tttrace -trace of the tensor train operator

## Implementation structure

- Computes the trace of a tensor train operator. Syntax:
- tttrace=trace(tt)
- tt -tensor train operator
- tttrace -trace of the tensor train operator
- Read sizes and ranks
- Make an auxiliary tensor train
- Run through all tensor trains
- Preallocate a core
- Fill in the core
- Reshape the core
- Sum up the auxiliary tensor train
- Pronouncement of experts to the effect that something
