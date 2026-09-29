# kernel/utilities/killcross.m

## Purpose

`killcross.m` zeroes the specified rows and columns of a matrix. Syntax:

```
M=killcross(M,f1idx,f2idx)
```

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/killcross.m>

## Behaviour

- The function calls the internal consistency checker `grumble(M,f1idx,f2idx)` before modifying the matrix.
- It then wipes the specified indices: `M(f2idx,:)=0; M(:,f1idx)=0;`, i.e. all rows listed in `f2idx` and all columns listed in `f1idx` are set to zero.
- Consistency enforcement (`grumble`) raises errors when:
  - `M` is not numeric or not a matrix (`'M must be a matrix.'`).
  - `f1idx` or `f2idx` is not numeric, not real, contains values less than 1, or contains non-integer values (`'index arrays must contain positive integers.'`).
  - `f1idx` or `f2idx` contains repeated elements (`'repeated elements not allowed in the index arrays.'`).
  - Any element of `f1idx` exceeds `size(M,2)` or any element of `f2idx` exceeds `size(M,1)` (`'index array element exceeds spectrum matrix dimension.'`).

## Inputs and outputs

Inputs:

- `M` — a matrix.
- `f1idx` — numbers of the columns that should be zeroed.
- `f2idx` — numbers of the rows that should be zeroed.

Outputs:

- `M` — a matrix with the specified rows and columns zeroed.

## References

- Spinach Wiki: <https://spindynamics.org/wiki/index.php?title=killcross.m>
- Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/killcross.m>
