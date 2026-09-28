# kernel/utilities/dipolar.m

- Signature: `spin_system=dipolar(spin_system)`

## Purpose

Computes dipolar couplings with or without periodic boundary conditions. This is an auxiliary Spinach kernel function; direct calls are discouraged. Use `xyz2dd` and `xyz2hfc` to convert Cartesian coordinates into dipolar and hyperfine couplings, respectively.

## Physical / mathematical content

For each qualifying distance vector, the function computes a dipolar coupling matrix from the two spins’ gyromagnetic ratios, their separation, and the unit vector along that separation. The prefactor includes a factor of `0.5` to account for double counting.

## Numerical / algorithmic content

- Examines distinct, coordinate-specified spin pairs within each chemical subsystem. With periodic boundaries, it considers images across one, two, or three translation directions, using the range set by `spin_system.tols.dd_ncells`.
- Retains distance vectors shorter than `spin_system.tols.prox_cutoff` and raises an error for separations below `0.5` Angstrom. It records qualifying pairs in `spin_system.inter.proxmatrix`.
- Optionally applies the `sodd` spin-orbit correction, removes the matrix trace to clean up numerical noise, and adds each contribution to the corresponding coupling matrix.

## Parameters / inputs

- spin_system -Spinach data object containing infor-
- mation about chemical subsystems, ato-
- mic coordinates, and periodic bounda-
- ry conditions

## Outputs

- spin_system -Spinach data object with the interac-
- tion arrays updated with dipolar and
- hyperfine coupling information

## Implementation structure

The function checks that the spin-system object has essential fields, reports the distance threshold and number of qualifying spin pairs, then accumulates dipolar coupling matrices for their retained distance vectors.