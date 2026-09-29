# kernel/overloads/@opium/mtimes.m

- Signature: `c=mtimes(a,b)`

## Meaning and behaviour

An `opium` object represents the scaled unit matrix `coeff*I_dim`. The overload keeps that representation for scalar scaling and for products of two `opium` objects; it returns a numeric result for the nonscalar numeric branches.

- Numeric scalar `a` times `opium` object `b`: the first branch copies `b` and replaces its `coeff` by `a*b.coeff`. The branch checks that `a` is a numeric scalar; it does not separately test `b` with `isa(...,'opium')` before accessing `b.coeff`.
- `opium` object `a` times numeric scalar `b`: copies `a` and replaces its coefficient by `b*a.coeff`; the scalar remains implicit in the object.
- Numeric nonscalar `a` times `opium` object `b`: checks `size(a,2)==size(b,1)`, then returns `b.coeff*a`.
- `opium` object `a` times numeric nonscalar `b`: checks `size(a,2)==size(b,1)`, then returns `a.coeff*b`.
- Two `opium` operands: checks the same size equality, then copies `a` and sets its coefficient to `a.coeff*b.coeff`. Since each represented matrix is square, the check requires their reported dimensions to match; the result retains `a`'s dimension.

In the nonscalar numeric cases, the represented scaled identity acts as multiplication by its coefficient after the stated conformability check; these branches do not construct the identity with `eye` or `speye`. The source checks only the displayed pair of dimensions in these branches and supplies no cell-array or broadcasting rule. Other inputs that reach the final branch raise an error; the source's predicates and field accesses, rather than a broader inferred type contract, define the branches above.

## Inputs and output

- Inputs: numeric scalar/matrix operands and/or `opium` objects on the branches described above.
- Output: a coefficient-updated `opium` for scalar scaling or an `opium`-by-`opium` product; a numeric coefficient-scaled operand for the nonscalar numeric branches.

## Source links

- MATLAB source: [kernel/overloads/@opium/mtimes.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/overloads/@opium/mtimes.m)
- Existing Wiki page: [opium/mtimes.m](https://spindynamics.org/wiki/index.php?title=opium/mtimes.m)
- DOI retained from a source-file comment: [10.1063/1.1719961](https://doi.org/10.1063/1.1719961)
