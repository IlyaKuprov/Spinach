# kernel/overloads/@ttclass/ttclass.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/ttclass.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=ttclass.m)

## Signature

`tt=ttclass(coeff,kronterms,tolerance)`

## Purpose and representation

Constructs a tensor-train object from Kronecker-product factors. Each column of the two-dimensional `kronterms` cell array describes one train: row `d` is factor/core `d`, so the cell array is ordered as cores by rows and buffered trains by columns. Each factor matrix is converted to full storage and reshaped to `[1,size(factor,1),size(factor,2),1]`, giving singleton left and right boundary ranks. Multiple columns represent the sum of the corresponding trains, with one coefficient and one tolerance per column.

The tensor-train format stores high-dimensional Kronecker products compactly; see [the cited tensor-train reference](https://doi.org/10.1137/090752286). The constructor stores the supplied tolerance as metadata; it does not itself truncate or recompress the factors.

## Inputs

- `coeff` — numeric row vector of coefficients; complex values are permitted.
- `kronterms` — nonempty two-dimensional cell array whose entries are matrices. Its column count must equal `numel(coeff)`.
- `tolerance` — real, non-negative numeric row vector, with one entry per coefficient. The source header describes each entry as the maximum allowed 2-norm deviation between the tensor-train and flat-matrix representations.

Call with no inputs to receive the class's default-initialised object; otherwise the constructor requires all three inputs. It sets `tt.debuglevel=0`.

## Output and indexing

- `tt` — a `ttclass` object with `tt.coeff=coeff`, `tt.tolerance=tolerance`, and `tt.cores` shaped as `number of cores` by `number of buffered trains`.
- `tt.ncores` is the number of rows in `tt.cores`; `tt.ntrains` is its number of columns.
