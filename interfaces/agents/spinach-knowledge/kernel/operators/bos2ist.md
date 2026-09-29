# kernel/operators/bos2ist.m

- Signature: `[states,coeffs]=bos2ist(prod_spec,nlevels)`

## Meaning

Expands an ordered product of truncated-mode bosonic operators into Spinach irreducible spherical-tensor (IST) basis states and their coefficients. The accepted symbols in `prod_spec` are `C` (creation), `A` (annihilation), and `N` (number operator); `CCAA` is the source's example. Starting with the sparse identity, the routine scans the string from left to right and right-multiplies by the corresponding `weyl(nlevels)` matrix at each character. An empty character string leaves the identity matrix to be expanded.

## Index mapping and output

The routine passes the completed matrix to `oper2ist`. Its `states` are Spinach IST linear basis indices, not oscillator population labels; `lin2lm` converts an index to spherical-tensor `L,M` labels. The parallel `coeffs` values give the expansion coefficients for those states. In the called `oper2ist` implementation, linear labels start at zero and terms with coefficient magnitude no greater than `10*eps('double')` are omitted.

## Inputs and guards

`prod_spec` must be a character value whose characters are all in `C`, `A`, or `N`. The local check requires `nlevels` to be numeric, real, scalar, and at least one; the subsequent `weyl` call additionally enforces a positive integer. The truncation sets the matrix size used for the Weyl operators and IST expansion.

## References

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/operators/bos2ist.m)
- [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=bos2ist.m)
