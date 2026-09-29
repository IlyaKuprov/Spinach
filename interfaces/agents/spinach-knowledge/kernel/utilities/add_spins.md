# kernel/utilities/add_spins.m

## Purpose

Reduces the direct product of two su(2) irreducible representations into a direct sum of irreducible representations, returning the multiplicities of the total spin values that occur and the corresponding projection operators.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/add_spins.m>

## Behaviour

- Syntax: `[mult,proj]=add_spins(spin_a,spin_b)`.
- Validates both inputs through an internal consistency check (`grumble`): each must be numeric, real, scalar, at least `1/2`, and such that `2*spin+1` is an integer; otherwise an error is thrown (`'spin_a must be a positive integer or half-integer.'` / `'spin_b must be a positive integer or half-integer.'`).
- Builds the individual spin irreps via `pauli(2*spin+1)` for each input.
- Constructs the direct-product representation generators:
  - `Sx=kron(spin_a.u,spin_b.x)+kron(spin_a.x,spin_b.u)`
  - `Sy=kron(spin_a.u,spin_b.y)+kron(spin_a.y,spin_b.u)`
  - `Sz=kron(spin_a.u,spin_b.z)+kron(spin_a.z,spin_b.u)`
- Diagonalises the Casimir operator `Sx^2+Sy^2+Sz^2` and indexes its eigenvalues with `unique(uint32(D))` to identify the distinct total-spin sectors.
- Forms one projector per distinct eigenvalue from the eigenvector blocks of the Casimir diagonalisation.
- Canonicalises each projector block:
  - Records the multiplicity `mult(n)` as the number of columns of the projector.
  - Diagonalises the projected `Sz` block, sorts eigenvalues in descending order, and rotates the projector so that `Sz` is diagonal; fails with `'irrep canonicalisation failed.'` if the projected `Sz` does not match the canonical `pauli` `z` matrix within `sqrt(eps)` in the 1-norm.
  - Applies column sign flips until the projected `Sx` block has real positive entries (within `sqrt(eps)`); fails with the same error if the projected `Sx` or `Sy` blocks do not match the canonical `pauli` `x`/`y` matrices within `sqrt(eps)` in the 1-norm.

## Inputs and outputs

**Inputs**

- `spin_a` — quantum number of the first spin; an integer or a half-integer.
- `spin_b` — quantum number of the second spin; an integer or a half-integer.

**Outputs**

- `mult` — multiplicities corresponding to the values of the total spin that are present.
- `proj` — projectors that reduce the direct product representation; a cell array of matrices.

## References

- Spinach Wiki: <https://spindynamics.org/wiki/index.php?title=add_spins.m>
- Source file: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/add_spins.m>
