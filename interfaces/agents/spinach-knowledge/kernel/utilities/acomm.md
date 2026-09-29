# kernel/utilities/acomm.m

## Purpose

`acomm.m` computes the anticommutator of two square matrices, returning `C = A*B + B*A`. It is a simple shorthand provided by the Spinach kernel utilities.

Source: [kernel/utilities/acomm.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/acomm.m)

## Behaviour

- The function enforces input consistency before computing the result:
  - `A` must be a numeric square matrix, otherwise it errors with `'A must be a numeric square matrix.'`.
  - `B` must be a numeric square matrix, otherwise it errors with `'B must be a numeric square matrix.'`.
  - `A` and `B` must have identical dimensions, otherwise it errors with `'A and B must have the same dimensions.'`.
- Checks use `isnumeric`, `ismatrix`, and comparison of `size(A,1)` with `size(A,2)` (and likewise for `B`), plus `isequal(size(A),size(B))`.
- After validation, the function returns `C = A*B + B*A`.

## Inputs and outputs

**Syntax**

```matlab
C = acomm(A,B)
```

**Inputs**

- `A` — a numeric square matrix.
- `B` — a numeric square matrix of the same dimensions as `A`.

**Outputs**

- `C` — a square matrix, the anticommutator `A*B + B*A`.

## References

- Spin Dynamics Wiki page for this function: <https://spindynamics.org/wiki/index.php?title=acomm.m>
- Source file: [kernel/utilities/acomm.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/acomm.m)
