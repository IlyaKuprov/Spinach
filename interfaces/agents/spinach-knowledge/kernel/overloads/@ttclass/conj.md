# kernel/overloads/@ttclass/conj.m

## Links

- [Direct MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/conj.m)
- [Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=ttclass/conj.m)

## Storage and scope

In this `ttclass` storage, `tt.cores` is an `ncores`-by-`ntrains` cell array. Core `tt.cores{k,n}` has left/right bond-rank axes 1 and 4 and physical row/column axes 2 and 3. Adjacent cores contract by summing over their matching right/left bond index; each train has boundary ranks one. The row coefficient `tt.coeff(1,n)` weights train `n`, and the columns store separate coefficient-weighted TT chains. This is a tensor-train, not a `polyadic`, representation.

## Signature

`tt=conj(tt)`

## Behaviour

The function loops over every train and core and applies complex conjugation elementwise to each core, then applies complex conjugation to the coefficient array. It changes no core ordering, bond ranks, physical mode sizes, train count, or matrix/vector dimensions; the output has the same tensor-train layout.

## Checks

The overload contains no input, shape, rank, or coefficient validation and no error branch; it assumes a valid `ttclass` object.
