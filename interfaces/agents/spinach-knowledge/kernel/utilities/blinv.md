# kernel/utilities/blinv.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/blinv.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/blinv.m)

## Purpose

Computes Blicharski's relaxation theory invariants for a real 3x3 interaction matrix, as given by Equations 20-21 in [http://doi.org/10.1515/zna-1972-1012](http://doi.org/10.1515/zna-1972-1012).

## Behaviour

- Validates the input via an internal consistency check (`grumble`), which errors with `'A must be a real 3x3 matrix.'` if the argument is not numeric, not real, not a matrix, or not of size 3x3.
- Computes the first rank invariant `Lsq` as the sum of squared antisymmetric parts:
  - `(A(1,2)-A(2,1))^2 + (A(1,3)-A(3,1))^2 + (A(2,3)-A(3,2))^2`
- Computes the second rank invariant `Dsq` as:
  - `A(1,1)^2 + A(2,2)^2 + A(3,3)^2 - A(1,1)*A(2,2) - A(1,1)*A(3,3) - A(2,2)*A(3,3) + (3/4)*((A(1,2)+A(2,1))^2 + (A(1,3)+A(3,1))^2 + (A(2,3)+A(3,2))^2)`
- The function is not sensitive to the trace of the matrix.

## Inputs and outputs

**Inputs:**

- `A` — a real 3x3 matrix (the interaction matrix).

**Outputs:**

- `Lsq` — first rank invariant.
- `Dsq` — second rank invariant.

## References

- Blicharski, J. S. — Equations 20-21, [http://doi.org/10.1515/zna-1972-1012](http://doi.org/10.1515/zna-1972-1012).
- Spinach Wiki: [https://spindynamics.org/wiki/index.php?title=blinv.m](https://spindynamics.org/wiki/index.php?title=blinv.m)
