# kernel/overloads/@ttclass/diag.m

## Links

- [Direct MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/diag.m)
- [Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=ttclass/diag.m)

## Storage and scope

In this `ttclass` storage, `tt.cores` is an `ncores`-by-`ntrains` cell array. Core `tt.cores{k,n}` has left/right bond-rank axes 1 and 4 and physical row/column axes 2 and 3. Adjacent cores contract by summing over their matching right/left bond index; each train has boundary ranks one. The row coefficient `tt.coeff(1,n)` weights train `n`, and the columns store separate coefficient-weighted TT chains. This is a tensor-train, not a `polyadic`, representation.

## Signature

`tt=diag(tt)`

## Behaviour and result dimensions

The source obtains mode sizes with `sizes(tt)` and bond ranks with `ranks(tt)`. A vector is selected when all row-mode sizes are one or all column-mode sizes are one. For each train and core, it forms a square physical core by applying MATLAB `diag` to the reshaped physical slice for every left/right bond-index pair. The output represents a diagonal matrix whose side length is the product of the vector mode sizes.

Otherwise, the source accepts a matrix only when every mode is square (`sz(:,1)==sz(:,2)`). It extracts each core's physical diagonal separately for every bond-index pair and stores a vector core with physical dimensions `d-by-1`. The output is a column vector with length equal to the product of the input mode sizes. Rank bonds and train count are retained in both branches. Any input that is neither a recognised vector nor modewise square raises `Input should be either a square matrix or a vector.`

## Checks and limits

The vector/matrix decision is based on the reported mode sizes; `sizes(tt)` obtains these sizes from the first train, while the core transformations loop over all trains. This overload has no separate per-train mode-size or bond-compatibility check.
