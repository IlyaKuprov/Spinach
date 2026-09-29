# kernel/overloads/@polyadic/prefix.m

[Source on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@polyadic/prefix.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=polyadic/prefix.m)

## Signature

`p=prefix(a,p)`

## Operation and order

This overload adds a matrix on the left of the represented polyadic, so the intended product is `a*P`. If `a` is nonscalar, it checks `size(a,2)==size(p,1)` and prepends `a` to `p.prefix` (`p.prefix=[{a} p.prefix]`). Repeated calls therefore store the newest left factor first: a later `prefix(b,p)` represents `b*(a*P)`. The factor is recorded in the composition array rather than multiplying through the Kronecker terms. The source comments allow `a` itself to be polyadic.

If `a` is scalar, the method instead multiplies it directly into the first factor of every term, `p.cores{n}{1}=a*p.cores{n}{1}`; it does not append the scalar to `p.prefix`. This is a representation-level update using the ordinary overloaded `*` operation.

## Inputs, checks, and output

- `p` must satisfy `isa(p,'polyadic')`, or the method raises `p must be polyadic.`
- For nonscalar `a`, its column count must equal `size(p,1)`; otherwise the method raises `matrix dimension mismatch.` There is no separate input-type check for `a` in this function.
- `p`: the updated polyadic object. A nonscalar left factor sets the row count from `size(a,1)` while leaving the represented column dimension unchanged.

The function does not conjugate or transpose `a`, and it defines no broadcasting rule; only the scalar branch and the nonscalar dimension check above are explicit here. Related pages: [polyadic representation](./polyadic.md) and [size](./size.md).