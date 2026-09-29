# kernel/utilities/add_spins.m

## Purpose

Reduces the direct product of two su(2) irreducible representations into a direct sum of irreducible representations, returning the multiplicities of the total spin values that occur and the corresponding projection operators.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/add_spins.m>

## How to use it

`[mult,proj]=add_spins(spin_a,spin_b)` reduces the tensor product of two spins into total-spin sectors. Each input must be a real scalar integer or half-integer quantum number of at least 1/2; invalid inputs raise an error. The sectors are ordered by increasing total-spin Casimir eigenvalue. `mult` records the dimension of each sector, and `proj{n}` gives its basis vectors as columns in the original direct-product space.

The projected spin generators are canonicalised to the standard spin matrices: within each sector, columns follow descending `Sz` eigenvalue and phase conventions chosen to match `Sx` and `Sy`. The function errors if that canonicalisation fails, rather than returning an inconsistent projector. For spin quantum numbers a and b, the expected sectors have total spin from |a−b| to a+b in integer steps.

## Inputs and outputs

**Inputs**

- `spin_a` — quantum number of the first spin; an integer or a half-integer.
- `spin_b` — quantum number of the second spin; an integer or a half-integer.

**Outputs**

- `mult` — one value per total-spin sector; the implementation records the dimension of each projected block.
- `proj` — projectors that reduce the direct product representation; a cell array of matrices.

## References

- Spinach Wiki: <https://spindynamics.org/wiki/index.php?title=add_spins.m>
- Source file: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/add_spins.m>
