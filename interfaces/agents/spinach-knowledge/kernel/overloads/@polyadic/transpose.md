# kernel/overloads/@polyadic/transpose.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@polyadic/transpose.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=polyadic/transpose.m)

- Signature: `p=transpose(p)`

## Purpose

Returns the ordinary, non-conjugating transpose of the matrix represented by a polyadic, transforming its stored factorisation rather than materialising the full matrix.

## Factor and operator order

For every buffered term, the method applies MATLAB `transpose` to each stored core matrix in place in the nested cell representation; it does not reverse the core-factor cell order. It reverses the boundary chains while transposing their members: the new prefix is the reversed old suffix with each factor transposed, and the new suffix is the reversed old prefix with each factor transposed. This implements the product rule `(AB)^T = B^T A^T` for those ordered boundary products. The represented matrix dimensions are therefore exchanged.

## Inputs and output

- The signature accepts `p` and returns the transformed `p`.
- This implementation has no explicit class or consistency check; it directly accesses `p.cores`, `p.suffix`, and `p.prefix`.
- It uses ordinary transpose, not conjugate transpose, so complex entries are not conjugated.

## Materialisation

Each stored core and boundary factor is transformed eagerly, but the factorised representation is retained; no full matrix is formed.
