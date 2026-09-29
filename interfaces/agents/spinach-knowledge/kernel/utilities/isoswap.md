# kernel/utilities/isoswap.m

## Purpose

Makes isotope replacements in the input structures. All interactions are automatically scaled as necessary.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/isoswap.m>

## Behaviour

- Syntax: `[sys,inter]=isoswap(sys,inter,spins,new_iso)`.
- Enforces input consistency via an internal `grumble` subfunction, which errors if isotope information is missing from `sys`, if `inter` is not a structure, if `spins` is not a real numeric vector of integer indices not exceeding the number of isotopes, or if `new_iso` is not a character string.
- Wipes quadratic couplings before replacement, in both eigensystem representation (`inter.coupling.eigs`, together with `inter.coupling.euler`) and matrix representation (`inter.coupling.matrix`), printing a warning per affected spin that the coupling is not transferable.
- Wipes high-rank couplings stored in `inter.giant.coeff` and `inter.giant.euler`, printing a warning per affected spin.
- For each specified spin, computes `gamma_ratio = spin(new_iso)/spin(sys.isotopes{n})` and scales all couplings to every other spin `k` (both `{n,k}` and `{k,n}` entries) in the eigensystem and matrix representations by this ratio.
- Replaces the isotope string `sys.isotopes{n}` with `new_iso` for each specified spin.
- Quadratic and higher order couplings are wiped with a warning because they are not transferable.

## Inputs and outputs

Inputs:

- `sys`, `inter` — Spinach input data structures.
- `spins` — a vector of integers specifying spin numbers.
- `new_iso` — character string specifying the new isotope.

Outputs:

- `sys`, `inter` — Spinach input data structures with the isotope replacement applied and interactions scaled.

## References

- Spinach Wiki: <https://spindynamics.org/wiki/index.php?title=isoswap.m>
- Source file: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/isoswap.m>
