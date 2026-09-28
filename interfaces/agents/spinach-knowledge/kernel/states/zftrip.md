# kernel/states/zftrip.m

- Signature: `rho=zftrip(spin_system,ZFS,pops,Z,B,idx)`

## Purpose

Projects a zero-field triplet state with user-specified Cartesian ZFS-eigenstate populations onto the eigenstates of the higher-field ZFS plus Zeeman Hamiltonian. This is used, for example, for photo-generated two-electron triplets in triplet DNP. Syntax: `rho=zftrip(spin_system,ZFS,pops,Z,B,idx)`.

## Physical / mathematical content

Constructs the spin-1 ZFS Hamiltonian from the laboratory-frame tensor `ZFS`, diagonalizes it, and assigns the `pops` entries to its `X`, `Y`, and `Z` eigenstates. It then combines the ZFS and field-dependent Zeeman Hamiltonians, removes coherences in that high-field eigenbasis, and converts the resulting density matrix to Spinach spherical-tensor states for triplet `idx`.

## Numerical / algorithmic content

Two Hermitian eigenproblems are solved: first for the zero-field ZFS Hamiltonian, then for the ZFS plus Zeeman Hamiltonian at field `B`. The zero-field eigenvectors are ordered by increasing absolute eigenvalue to apply the organic triplet convention.

## Parameters / inputs

- `ZFS` - 3x3 ZFS tensor (Hz) in the laboratory frame of reference; use `zfs2mat()` to get it from D, E, and molecular Euler angles.
- `pops` - three-element vector with populations of the X, Y, and Z eigenstates of the ZFS tensor at zero magnetic field, order: `[pX pY pZ]`. X, Y, and Z are labelled using the organic triplet convention `|Dzz|>|Dxx|>|Dyy|`, under which D and E have opposite signs, `-1/3<E/D<0` (Poole, Farach, Jackson, J. Chem. Phys. 61, 2220 (1974), DOI 10.1063/1.1682294); populations quoted in the transition metal convention `|Dzz|>|Dyy|>|Dxx|` with `0<E/D<1/3` must have X and Y swapped before the call; at E=0 the X and Y states are degenerate (all three when D=0 as well) and the labelling within the degenerate set is undefined; at E/D=-1/3 the X and Z energies are opposite in sign, not degenerate, but equal in magnitude, so the sort by |energy| cannot tell X from Z; in both cases the affected populations must be equal for the result to be meaningful.
- `Z` - 3x3 Zeeman interaction tensor (Hz/Tesla) in the laboratory frame of reference; use functions like `axrh2mat()` to get it from eigenvalues and molecular Euler angles.
- `B` - magnetic field directed along the Z axis of the laboratory frame of reference, Tesla.
- `idx` - index of the electron triplet (use `'E3'`) in the `sys.isotopes` list.

## Outputs

- `rho` - spin density matrix (Hilbert space) or state vector (Liouville space).

## Implementation structure

- Checks tensor symmetry, population normalization, field, and triplet-spin index; builds the zero-field density matrix from the labelled populations, diagonalizes the high-field Hamiltonian, discards high-field coherences, and expands the retained density operator in spherical-tensor basis states.
