# kernel/overloads/@ttclass/clearcoeff.m

## Links

- [Direct MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/clearcoeff.m)
- [Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=ttclass/clearcoeff.m)

## Storage and scope

In this `ttclass` storage, `tt.cores` is an `ncores`-by-`ntrains` cell array. Core `tt.cores{k,n}` has left/right bond-rank axes 1 and 4 and physical row/column axes 2 and 3. Adjacent cores contract by summing over their matching right/left bond index; each train has boundary ranks one. The row coefficient `tt.coeff(1,n)` weights train `n`, and the columns store separate coefficient-weighted TT chains. This is a tensor-train, not a `polyadic`, representation.

## Signature

`tt=clearcoeff(tt)`

## Behaviour

For each train `n`, the function computes `tt.coeff(1,n)^(1/ncores)`, multiplies every core in that train by this factor, then sets that coefficient to one. Applying the same factor to all `ncores` cores absorbs the train coefficient into the core chain; the represented value is unchanged. Core count, train count, mode sizes, ranks, and output shape are unchanged.

## Checks

The overload contains no input, shape, rank, or coefficient validation and no error branch; it assumes a valid `ttclass` object with consistent core and coefficient storage.
