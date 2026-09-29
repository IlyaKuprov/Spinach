# kernel/utilities/remncomm.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/remncomm.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/remncomm.m)

## Purpose

Removes from a Hermitian operator `A` the part that does not commute with a Hermitian operator `B`, returning the commuting part `C`. The operator `B` is supplied through its eigenvectors and eigenvalues rather than as a matrix.

## Behaviour

- Syntax: `C=remncomm(A,EvecB,EvalB)`.
- The function first validates its inputs via an internal consistency check (`grumble`).
- `A` is transformed into the eigenbasis of `B` as `EvecB'*A*EvecB`.
- In that basis, elements of `A` that link eigenvalues of `B` differing by more than eigensolver roundoff are zeroed out. The degeneracy mask is `abs(EvalB-EvalB.')<=100*numel(EvalB)*eps(max(EvalB)-min(EvalB))`, and the masked matrix is `A.*degen_mask`.
- The commuting part is transformed back to the original basis as `EvecB*A*EvecB'`.
- Within a degenerate eigenspace of `B`, every Hermitian operator supported on that eigenspace commutes with `B`, so the corresponding block of `A` (not just its diagonal) is kept.

Input validation errors:

- `A` must be a Hermitian matrix (numeric, square, and `ishermitian`).
- `EvecB` must be a square numeric array of column vectors.
- `EvalB` must be a finite real floating-point column vector with as many elements as `EvecB` has columns.
- The spread `max(EvalB)-min(EvalB)` must be representable in floating point.

## Inputs and outputs

**Inputs**

- `A` — a square (Hermitian) matrix.
- `EvecB` — a square matrix containing the eigenvectors of `B` in columns.
- `EvalB` — a column vector containing the eigenvalues of `B` in the same order as the columns of `EvecB`.

**Outputs**

- `C` — a square matrix: the part of `A` that commutes with `B`.

## References

- Spinach Wiki: [remncomm.m](https://spindynamics.org/wiki/index.php?title=remncomm.m)
