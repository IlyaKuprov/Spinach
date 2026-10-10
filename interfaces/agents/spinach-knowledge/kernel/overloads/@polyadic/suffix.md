# kernel/overloads/@polyadic/suffix.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@polyadic/suffix.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=polyadic/suffix.m)

- Signature: `p=suffix(p,a)`

## Purpose

Adds a right-boundary factor to a factorised polyadic representation. The operation edits its stored factors; it does not expand the represented matrix.

## Representation and operator order

The polyadic stores buffered terms in the outer cell array `p.cores`; each term is a cell sequence of core matrices. Boundary factors are stored in ordered cell sequences `p.prefix` and `p.suffix`. The source documentation notes that a suffix may itself be a polyadic.

For non-scalar `a`, the method appends `a` to `p.suffix`, so it remains a separate factor rather than being multiplied into the cores. A suffix acts on an argument before the polyadic does; appending a factor extends the right-boundary product chain. The compatibility check is `size(p,2) == size(a,1)`.

A numeric scalar `a` is represented as `opium(size(p,2),a)` and appended through the same suffix path. This dimensioned scaled identity preserves both matrix dimensions and the polyadic return type, including for zero. A scalar-sized polyadic remains a matrix factor, not a numeric scaling coefficient.

## Inputs and output

- `p` must be a `polyadic` object; the method checks this before doing either branch.
- Numeric scalar coefficients are distinguished from polyadic factors. A row-dimension mismatch raises `matrix dimension mismatch.`; otherwise the factor is appended.
- Returns the updated `p`. The method does not add a separate explicit check of `a`'s type.

## Materialisation

All factors remain in the suffix chain, so this overload does not form a dense polyadic matrix. Explicit `simplify` may absorb scaled identities through scalar multiplication or return a numeric zero matrix.
