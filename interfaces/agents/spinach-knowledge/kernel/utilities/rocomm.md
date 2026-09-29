# kernel/utilities/rocomm.m

## Purpose

`rocomm.m` computes the right-ordered nested commutator `[[[[A{1},A{2}],A{3}],A{4}],...]` from a user-supplied cell array of matrices. It is part of the Spinach kernel utilities.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/rocomm.m>

## Behaviour

- Syntax: `C=rocomm(A)`.
- The function first validates the input via an internal consistency check (`grumble`), which errors with `'A must be a cell array of square matrices.'` if `A` is not a cell array, if any element is not numeric, or if any element is not square (row count differs from column count).
- The nesting starts with `C=A{1}` and iterates `n=2:numel(A)`, updating `C=C*A{n}-A{n}*C` at each step. This builds the right-ordered nested commutator by successively commuting the running result with each subsequent matrix in the cell array.

## Inputs and outputs

**Inputs**

- `A` — a cell array of square matrices.

**Outputs**

- `C` — the right-ordered nested commutator.

## References

- Spinach Dynamics Wiki page for `rocomm.m`: <https://spindynamics.org/wiki/index.php?title=rocomm.m>
