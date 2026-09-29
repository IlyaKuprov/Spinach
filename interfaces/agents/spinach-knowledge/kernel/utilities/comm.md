# kernel/utilities/comm.m

## Purpose

`comm.m` computes the commutator of two square matrices, providing a simple shorthand for the expression `A*B-B*A`. Source: [kernel/utilities/comm.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/comm.m).

## Behaviour

The function is called as `C=comm(A,B)`. It first runs an internal consistency check (`grumble`) on the inputs and then returns `C=A*B-B*A`. The consistency check raises an error with the message `both inputs must be numeric.` if either input is not numeric, and an error with the message `both inputs must be square matrices.` if either input is not square (i.e., its number of rows differs from its number of columns).

## Inputs and outputs

Inputs:
- `A`, `B` — square matrices; both must be numeric.

Outputs:
- `C` — a square matrix, the commutator `A*B-B*A`.

## References

- Spinach Wiki: <https://spindynamics.org/wiki/index.php?title=comm.m>
- Source file: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/comm.m>
