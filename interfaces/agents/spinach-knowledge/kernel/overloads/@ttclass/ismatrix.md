# kernel/overloads/@ttclass/ismatrix.m

- Signature: `answer=ismatrix(tt)`

## Purpose and behaviour

Returns logical true exactly when `tt` is a `ttclass` object and its `cores` property is non-empty; otherwise it returns logical false. The predicate checks the stored TT core-cell container only; it is not a test or description of a polyadic-object representation. It does not inspect core contents, validate TT ranks, or determine whether a general MATLAB array has exactly two dimensions, so read the result as this overload's storage predicate rather than a general shape validation.

## Sources

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@ttclass/ismatrix.m)
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=ttclass/ismatrix.m)
