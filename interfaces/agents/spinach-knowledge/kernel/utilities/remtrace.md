# kernel/utilities/remtrace.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/remtrace.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/remtrace.m)

## Purpose

Subtracts an appropriate multiple of the unit matrix from a square matrix to make it traceless.

## Behaviour

- Syntax: `A=remtrace(A)`.
- The function calls an internal consistency check (`grumble`) that errors with `'A must be a square matrix.'` if the input is not numeric or if it is not square (`size(A,1)~=size(A,2)`).
- The trace is removed by computing `dim=size(A,1)` and updating `A=A-speye(dim)*trace(A)/dim`, i.e. subtracting `trace(A)/dim` times the identity matrix of the same dimension.

## Inputs and outputs

**Input:**

- `A` — a square matrix.

**Output:**

- `A` — a square matrix with a zero trace.

## References

- Spinach Dynamics Wiki: [remtrace.m](https://spindynamics.org/wiki/index.php?title=remtrace.m)
