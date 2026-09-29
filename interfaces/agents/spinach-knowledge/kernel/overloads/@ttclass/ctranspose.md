# kernel/overloads/@ttclass/ctranspose.m

## Links

- [Direct MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/ctranspose.m)
- [Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=ttclass/ctranspose.m)

## Storage and scope

In this `ttclass` storage, `tt.cores` is an `ncores`-by-`ntrains` cell array. Core `tt.cores{k,n}` has left/right bond-rank axes 1 and 4 and physical row/column axes 2 and 3. Adjacent cores contract by summing over their matching right/left bond index; each train has boundary ranks one. The row coefficient `tt.coeff(1,n)` weights train `n`, and the columns store separate coefficient-weighted TT chains. This is a tensor-train, not a `polyadic`, representation.

## Signature

`ttrain=ctranspose(ttrain)`

## Behaviour

For every core, the function permutes dimensions with `[1 3 2 4]`, exchanging the row- and column-mode axes while preserving both bond-rank axes. It then calls `conj`, which conjugates every permuted core and all train coefficients. Thus an input matrix of mode-product dimensions `M-by-N` is represented as its Hermitian transpose, with dimensions `N-by-M`; ranks and core/train counts are unchanged.

## Checks

The overload performs no matrix-shape, rank, or class validation and has no explicit error branch. It applies the same dimension permutation to every stored core.
