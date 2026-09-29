# kernel/overloads/@polyadic/size.m

[Source on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@polyadic/size.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=polyadic/size.m)

## Signatures and results

- `size(p)` with zero or one requested output returns the two-element row vector `[nrows ncols]`.
- `[nrows,ncols]=size(p)` returns the row and column dimensions separately.
- `size(p,dim)` returns the requested dimension for `dim` equal to `1` (rows) or `2` (columns).

For rows, the method uses the row count of `p.prefix{1}` when a prefix exists; otherwise it multiplies the row dimensions of the factors in the first core term. For columns, it uses the column count of the last suffix factor when a suffix exists; otherwise it multiplies the column dimensions of the factors in the first core term. Thus the factorised dimensions are computed from stored factors without opening the Kronecker products. The method obtains dimensions from the first core term.

When `dim` is supplied it must be a scalar equal to `1` or `2`; otherwise the method raises `for a polyadic object, dim must be 1 or 2`. Call/output combinations that reach the function’s explicit fallback branch raise `invalid call syntax.` The source has no separate two-result branch when `dim` is supplied. This overload only inspects factor sizes; it does not multiply, conjugate, or broadcast factors. Related pages: [polyadic representation](./polyadic.md), [prefix](./prefix.md), and [simplify](./simplify.md).