# kernel/utilities/isoswap.m

## Purpose

Replaces selected isotopes in Spinach input structures. Use it after importing a calculation to obtain the bilinear couplings for a different isotope without changing the interpretation of the source data.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/isoswap.m>

## Behaviour

- Syntax: `[sys,inter]=isoswap(sys,inter,spins,new_iso)`.
- A selected source isotope must have nonzero magnetogyric ratio: its couplings cannot otherwise be rescaled by a ratio of magnetogyric ratios. A system without `inter.coupling` can still have its isotope labels replaced.
- Wipes quadratic couplings before replacement, in both eigensystem representation (`inter.coupling.eigs`, together with `inter.coupling.euler`) and matrix representation (`inter.coupling.matrix`), printing a warning per affected spin that the coupling is not transferable.
- Wipes high-rank couplings stored in `inter.giant.coeff` and `inter.giant.euler`, printing a warning per affected spin.
- For each selected spin, multiplies its bilinear couplings to every other spin by `spin(new_iso)/spin(old_iso)` in both pair directions. This applies to scalar couplings (`inter.coupling.scalar`), tensor matrices, and tensor eigenvalues. A pair with both isotopes replaced acquires one ratio from each partner.
- Replaces the isotope string `sys.isotopes{n}` with `new_iso` for each specified spin.

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
