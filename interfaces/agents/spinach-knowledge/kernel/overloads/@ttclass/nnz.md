# kernel/overloads/@ttclass/nnz.m

- Signature: `answer=nnz(ttrain)`
- Source: [`kernel/overloads/@ttclass/nnz.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/nnz.m)
- Wiki: [`ttclass/nnz.m`](https://spindynamics.org/wiki/index.php?title=ttclass/nnz.m)

## Core action

Applies MATLAB `nnz` to each entry of `ttrain.cores` and sums those counts. The scalar therefore counts stored nonzero core entries across the train buffer; it is not the number of nonzero entries in the represented matrix or tensor. It does not expand the train and does not use `ttrain.coeff`.

## Result and guards

Returns the summed count as a scalar. The method has no explicit type, shape, or rank guard, and performs no core or rank transformation. Its count follows the stored core shapes, not the logical dimensions represented by them.
