# kernel/utilities/isoswap.m

## Purpose

Replaces selected isotopes in Spinach input structures. Use it after importing a calculation to obtain the bilinear and the quadrupolar couplings for a different isotope without changing the interpretation of the source data.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/isoswap.m>

## Behaviour

- Syntax: `[sys,inter]=isoswap(sys,inter,spins,new_iso)`.
- A selected source isotope must have nonzero magnetogyric ratio: its couplings cannot otherwise be rescaled by a ratio of magnetogyric ratios. A system without `inter.coupling` can still have its isotope labels replaced.
- Rescales the quadratic coupling of each selected spin, in both eigensystem representation (`inter.coupling.eigs`, Euler angles kept) and matrix representation (`inter.coupling.matrix`), by the ratio of `Q/(2I(2I-1))` between the new and the old isotope, with quadrupole moments `Q` and spins `I` taken from `spin()`. Separate Gaussian 16 calculations for 35Cl/37Cl and 10B/11B are reproduced to better than 1e-4.
- When that ratio is not a finite nonzero number, the quadratic coupling is not transferable: it is wiped (together with `inter.coupling.euler` in eigensystem representation) and a warning is printed per affected spin. This happens when either isotope has no quadrupole moment (substituting 15N for 14N removes the quadrupolar coupling), when a moment is not tabulated, and for electron zero-field splittings.
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
