# kernel/utilities/hdot.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/hdot.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/hdot.m)

## Purpose

`hdot.m` computes the Frobenius inner product of two matrices via a Hadamard (element-wise) route, serving as an efficient replacement for `trace(A'*B)`.

## Behaviour

- The function evaluates `H = sum(conj(A).*B,'all')`, exploiting the identity `trace(A'*B) = sum(sum(conj(A).*B))`.
- This approach requires only O(n^2) multiplications, compared to O(n^3) for the direct `trace(A'*B)` computation.
- Before computing, a consistency check (`grumble`) enforces that both inputs are numeric and have identical dimensions, raising errors otherwise:
  - `'both inputs must be numeric.'` if either input is non-numeric.
  - `'the two inputs must have identical dimensions.'` if the sizes differ.

## Inputs and outputs

**Inputs:**

- `A`, `B` — square matrices of the same size; both must be numeric.

**Outputs:**

- `H` — the Frobenius inner product of `A` and `B`.

**Syntax:** `H=hdot(A,B)`

## References

- Spinach Wiki: [hdot.m](https://spindynamics.org/wiki/index.php?title=hdot.m)
