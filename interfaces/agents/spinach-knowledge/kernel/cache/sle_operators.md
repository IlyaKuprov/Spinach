# kernel/cache/sle_operators.m

- Signature: `[Lx,Ly,Lz,D,space_basis]=sle_operators(max_rank,int_ranks)`

## Purpose

Constructs the Wigner D function basis set and lab-space rotation generators required by the SLE module, along with product superoperators for requested interaction ranks.

## Physical / mathematical content

- `space_basis` indexes each basis Wigner function by `[L M N]`.
- `Lx`, `Ly`, and `Lz` represent lab-space rotation generators in this basis for building the lab-space diffusion operator.
- `D{r}{m,n}` represents multiplication by a rank-`r` Wigner function, for building the spin Hamiltonian operator.

## Numerical / algorithmic content

- Constructs the spatial basis for ranks `0:max_rank` and forms sparse rotation-generator matrices.
- Builds sparse product superoperators from Clebsch–Gordan coefficients and a Wigner-function normalization factor. Rank 2 uses hard-coded coefficient formulae; other ranks call `clebsch_gordan.m` and require the Java virtual machine during construction.
- Loads an existing operator set from a disk cache when available; otherwise attempts to save the constructed set. Cached sets load without the Java virtual machine.
- Uses `parfor` when computing Clebsch–Gordan coefficients.

## Parameters / inputs

- `max_rank` — maximum `L` rank for Wigner D functions; a positive integer.
- `int_ranks` — row vector of distinct positive interaction ranks for which product superoperators are required; may be empty if only rotation generators are needed.

## Outputs

- `space_basis` — lab-space basis descriptor in `[L M N]` format, giving the indices of each Wigner function.
- `Lx`, `Ly`, `Lz` — lab-space rotation generators in the Wigner function basis, used to build the lab-space diffusion operator.
- `D` — cell array indexed by interaction rank `r`; each requested `D{r}` is a `(2r+1)`-by-`(2r+1)` cell array of Wigner-function product superoperators used to build the spin Hamiltonian operator. `D{r}{m,n}` corresponds to multiplication by the Wigner function with projections `M=r+1-m` and `N=r+1-n`.

## Implementation structure

- Validates `max_rank` and `int_ranks`, then checks for a cache file keyed by both inputs.
- On a cache miss, enumerates the `[L M N]` basis, constructs `Lx`, `Ly`, and `Lz` from the raising generator, and builds product-action matrices for each requested interaction rank.
- Attempts to save the results to the cache; a write-protected installation produces a warning rather than preventing the results from being returned.

<https://spindynamics.org/wiki/index.php?title=sle_operators.m>
