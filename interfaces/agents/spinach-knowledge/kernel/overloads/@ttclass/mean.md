# kernel/overloads/@ttclass/mean.m

[Mapped MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/mean.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=ttclass/mean.m)

## Signature

`answer=mean(ttrain,dim)`

## Core and coefficient action

If every mode size is one, the method immediately returns `full(ttrain)` as a scalar, before choosing or validating `dim`. Otherwise, when `dim` is omitted it selects the first non-singleton dimension in matrix order (dimension 1 before dimension 2). Only dimensions 1 and 2 are handled; another value raises `incorrect dimension specificaton.`.

For every train and core, it sums the selected physical core dimension (array dimension `dim+1`), reshapes that dimension to one, and divides by that core's original size along the selected mode. For `dim=1`, each output core has shape `[rank(k,n),1,size2(k),rank(k+1,n)]`; for `dim=2`, it has shape `[rank(k,n),size1(k),1,rank(k+1,n)]`. Ranks and the other mode sizes are retained. The auxiliary train keeps `ttrain.coeff`; its tolerance is reset to `zeros(1,ntrains)`. There is no conjugation.

The result remains a `ttclass` unless all resulting mode sizes are one, in which case the method returns `full(answer)`. No full matrix is formed for a non-scalar result.

## Inputs and output

- `ttrain` — tensor-train representation of a matrix.
- `dim` — optional dimension, 1 or 2.
- `answer` — mean along the selected matrix dimension, as a tensor train or materialised scalar.
