# kernel/overloads/@polyadic/suffix.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@polyadic/suffix.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=polyadic/suffix.m)

- Signature: `p=suffix(p,a)`

## Purpose

Adds a right-boundary factor to a factorised polyadic representation. The operation edits its stored factors; it does not expand the represented matrix.

## Representation and operator order

The polyadic stores buffered terms in the outer cell array `p.cores`; each term is a cell sequence of core matrices. Boundary factors are stored in ordered cell sequences `p.prefix` and `p.suffix`. The source documentation notes that a suffix may itself be a polyadic.

For non-scalar `a`, the method appends `a` to `p.suffix`, so it remains a separate factor rather than being multiplied into the cores. A suffix acts on an argument before the polyadic does; appending a factor extends the right-boundary product chain. The compatibility check is `size(p,2) == size(a,1)`.

For scalar `a`, no suffix cell is added: the scalar eagerly left-multiplies the last core matrix in every buffered term, using `a * p.cores{n}{end}`. This branch has no separate dimension check.

## Inputs and output

- `p` must be a `polyadic` object; the method checks this before doing either branch.
- `a` is classified by `isscalar(a)`. For a non-scalar, a row-dimension mismatch raises `matrix dimension mismatch.`; otherwise the factor is appended.
- Returns the updated `p`. The method does not add a separate explicit check of `a`'s type.

## Materialisation

Scalar absorption eagerly scales each term's final stored core matrix. A non-scalar factor is retained in the suffix chain, so this overload does not form a dense polyadic matrix.
