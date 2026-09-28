# kernel/overloads/@ttclass/trace.m

- Signature: `tttrace=trace(tt)`

## Purpose

Compute the trace of a tensor-train operator by tracing each core's physical matrix dimensions and contracting the resulting train.

## Parameters / inputs

- `tt` — tensor train operator.

## Outputs

- `tttrace` — trace of the tensor-train operator.

## Implementation

The function reads the core sizes and ranks, then creates an auxiliary train with the same coefficients and zero tolerances. For every core and every pair of bond indices, it reshapes the physical dimensions to a matrix and stores that matrix's trace in a core with singleton physical dimensions. It converts the auxiliary train to its full value with `full` for the result.

## Source

D. Savostyanov and I. Kuprov, [`ttclass/trace.m`](https://spindynamics.org/wiki/index.php?title=ttclass/trace.m).
