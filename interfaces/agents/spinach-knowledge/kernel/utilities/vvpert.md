# kernel/utilities/vvpert.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/vvpert.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/vvpert.m)

## Purpose

Computes eigenvalue corrections and the Van Vleck transformation generator using Van Vleck perturbation theory, following Shavitt and Redmon, but excluding the quasi-degenerate split.

## Behaviour

- Syntax: `[Ep,G]=vvpert(E0,H1,order)`.
- Calls `grumble(E0,H1,order)` to enforce consistency of the inputs.
- Enforces exact Hermiticity of the perturbation via `H1=H1/2+H1'/2`.
- Builds reciprocal energy differences `Q=1./(E0'-E0)` with the diagonal zeroed using `speye`; if any element of `Q` is non-finite, it errors with `'H0 has degenerate energy levels.'`.
- Uses energy differences `delta=E0-E0.'` and Baker-Campbell-Hausdorff coefficients `bch=1./factorial(0:order)`.
- Allocates recursion storage `G`, `W` (each `order`-by-1 cell arrays) and `C` as an `(order+1)`-by-`(order+1)` cell array.
- For each perturbation order `n=1:order`:
  - Builds higher nested commutators independent of the current generator: for `m=2:n`, accumulates `K` over `k=(m-1):(n-1)` as `comm(C{m,k+1},G{n-k})` and stores it in `C{m+1,n+1}`.
  - Builds the transformed Hamiltonian term without the `H0` commutator: `K=H1` when `n==1`, otherwise `K=comm(H1,G{n-1})`; then adds `bch(m+1)*C{m+1,n+1}` for `m=2:n`.
  - Splits the term into the effective Hamiltonian `W{n}=diag(diag(K))` and the generator `G{n}=Q.*(K-W{n})`.
  - Stores the first nested commutator `C{2,n+1}=delta.*G{n}`, adding `comm(H1,G{n-1})` when `n>1`.
- Sums energy corrections: `Ep=E0` plus `real(diag(W{n}))` for each order `n`.
- Sums generator corrections by reshaping the `G` cell array into a 1-by-1-by-`numel(G)` array and summing over the third dimension with `cell2mat`.
- Notes from the header: no degeneracies are allowed in `H0`; `H1` must be finite and Hermitian; the theory only converges when `norm(H1,2)` is much smaller than the smallest energy gap in `H0`; complexity is cubic both in the order and in the matrix dimension; numerical artefacts appear beyond about order 10-12 for typical problems.

## Inputs and outputs

**Inputs**

- `E0` — eigenvalues of `H0`, a column vector of real numbers.
- `H1` — perturbation, written in the basis that diagonalises `H0`.
- `order` — order of perturbation theory to be used.

**Outputs**

- `Ep` — eigenvalues of `H0+H1` to the specified order, a column vector of reals, not necessarily sorted in the same way as the input.
- `G` — Van Vleck transformation generator, such that `expm(G)` is a square unitary matrix with eigenvectors in columns, in the same order as the eigenvalues in `Ep`.

**Validation (grumble)**

- Errors with `'E0 must be a real column vector.'` if `E0` is not numeric, not real, or not a column vector.
- Errors with `'H1 must be a finite Hermitian matrix.'` if `H1` is not numeric, not square, contains non-finite entries, or `norm(H1-H1',1)>1e-10*norm(H1,1)`.
- Errors with `'dimensions of E0 and H1 must be consistent.'` if `numel(E0)` does not match both dimensions of `H1`.
- Errors with `'order must be a positive real integer.'` if `order` is not numeric, not real, not scalar, not an integer, or less than 1.

## References

- Spinach Wiki: [vvpert.m](https://spindynamics.org/wiki/index.php?title=vvpert.m)
- Shavitt and Redmon (method followed by this implementation, excluding the quasi-degenerate split, as stated in the source header).
