# kernel/utilities/isoswap.m

- Signature: `[sys,inter]=isoswap(sys,inter,spins,new_iso)`

## Purpose

Replaces the isotope for each selected spin in the Spinach input structures and rescales transferable pairwise couplings.

## Parameters / inputs

- `sys` - Spinach system structure containing isotope specifications.
- `inter` - Spinach interaction structure.
- `spins` - vector of integer spin indices to replace.
- `new_iso` - character string specifying the replacement isotope.

## Outputs

- `sys`, `inter` - updated Spinach system and interaction structures.

## Numerical / algorithmic content

For each selected spin, the function computes the gyromagnetic-ratio factor `spin(new_iso)/spin(old_iso)`. It applies this factor to the spin's pairwise coupling entries in the eigensystem and/or matrix representations, when those fields are present, then updates the isotope string in `sys.isotopes`. Quadratic self-couplings in either representation and high-rank couplings are cleared because the source marks them as non-transferable; a warning is displayed when such terms are removed. The function checks that isotope data are present, `inter` is a structure, `spins` is a real numeric vector of integers with no entry greater than the isotope count, and `new_iso` is a character string.

## Source

[Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=isoswap.m)
