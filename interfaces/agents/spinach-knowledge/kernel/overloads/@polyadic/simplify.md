# kernel/overloads/@polyadic/simplify.m

[Source on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@polyadic/simplify.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=polyadic/simplify.m)

## Signature

`p=simplify(p)`

## Behaviour

The input must be a polyadic object; otherwise the function raises `p must be polyadic.` It first records `[nrows,ncols]=size(p)`. An empty core list, a zero prefix or suffix factor, or removal of all zero terms returns `spalloc(nrows,ncols,0)`, preserving the original represented dimensions as a sparse zero matrix.

Otherwise it repeatedly simplifies until a pass makes no changes. In `prefix` and `suffix`, nested polyadics are simplified recursively; zero factors collapse the result to the dimension-preserving sparse zero, identity factors are removed, and `opium` factors or numeric scalar factors are removed and their `coeff`/scalar is applied as `coeff*p` before returning.

Within each core term, nested polyadics are simplified recursively. A nested polyadic with no prefix or suffix and exactly one term is flattened by splicing its factor list into the parent term. A zero factor removes that whole summand. An identity factor that is not `opium` is replaced by `opium(size(factor,1),1)`. Empty terms are discarded; if none remain, the result is again the sparse zero of the saved dimensions. Adjacent `opium` factors are combined by an explicit `kron` call, and the process repeats.

Finally, if the result is still polyadic but has no prefix or suffix and contains exactly one term with exactly one factor, that factor is returned directly. These are local representation rewrites; the function does not expand the entire sum of Kronecker products. It does explicitly combine adjacent `opium` factors with `kron` and can return the sole factor instead of a polyadic wrapper.

No complex conjugation or broadcasting rule is implemented in this function. Related pages: [polyadic representation](./polyadic.md), [prefix](./prefix.md), and [size](./size.md).