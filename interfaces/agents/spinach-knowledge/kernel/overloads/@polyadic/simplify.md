# kernel/overloads/@polyadic/simplify.m

[Source on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@polyadic/simplify.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=polyadic/simplify.m)

## Signature

`p=simplify(p)`

## Behaviour

The input must be a polyadic object; otherwise the function raises `p must be polyadic.` It first records `[nrows,ncols]=size(p)`. An empty core list, a zero prefix or suffix factor, or removal of all zero terms returns `spalloc(nrows,ncols,0)`, preserving the original represented dimensions as a sparse zero matrix.

Otherwise it repeatedly simplifies until a pass makes no changes. In `prefix` and `suffix`, nested polyadics are simplified recursively; zero factors collapse the result to the dimension-preserving sparse zero, identity factors are removed, and `opium` factors or numeric scalar factors are removed and their `coeff`/scalar is applied as `coeff*p` before returning.

Within each core term, nested polyadics are simplified recursively. A nested polyadic with no prefix or suffix and exactly one term is flattened by splicing its factor list into the parent term. A zero factor removes that whole summand. An identity factor that is not `opium` is replaced by `opium(size(factor,1),1)`. Empty terms are discarded; if none remain, the result is again the sparse zero of the saved dimensions. Adjacent `opium` factors are combined by an explicit `kron` call, and the process repeats.

A singleton numeric core may be unwrapped. An implicit function-handle core retains its polyadic wrapper so that dimensions and adjoints remain available to norm estimation. Opaque cores are not tested for zero or identity. Flattening, term removal, and identity merging preserve their paired metadata.

No complex conjugation or broadcasting rule is implemented in this function. Related pages: [polyadic representation](./polyadic.md), [prefix](./prefix.md), and [size](./size.md).
