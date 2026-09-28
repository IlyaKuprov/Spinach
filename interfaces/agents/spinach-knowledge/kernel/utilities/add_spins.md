# kernel/utilities/add_spins.m

- Signature: `[mult,proj]=add_spins(spin_a,spin_b)`

## Purpose

Decompose the direct product of two `su(2)` spin representations into total-spin subspaces.

## Mathematical content

The function forms the combined `Sx`, `Sy`, and `Sz` generators as Kronecker sums, diagonalises the Casimir operator `Sx^2+Sy^2+Sz^2`, and groups its eigenvectors by eigenvalue. For each resulting subspace it orders the basis by descending eigenvalue of the projected z-spin operator and fixes basis-vector signs, then checks that the projected generators agree with the standard spin matrices from `pauli`.

## Parameters / inputs

- `spin_a` - positive integer or half-integer quantum number of the first spin.
- `spin_b` - positive integer or half-integer quantum number of the second spin.

## Outputs

- `mult` - number of basis vectors in each total-spin subspace found.
- `proj` - cell array whose matrices contain the corresponding subspace basis vectors, reducing the direct-product representation.