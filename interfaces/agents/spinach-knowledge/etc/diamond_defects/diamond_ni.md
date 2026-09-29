# etc/diamond_defects/diamond_ni.m

- MATLAB implementation: [etc/diamond_defects/diamond_ni.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/diamond_defects/diamond_ni.m)

- Signature: `[sys,inter]=diamond_ni(parameters)`
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=diamond_ni.m)
- W8 magnetic parameters: Isoya, Kanda, Norris, Tang, and Bowman, *Phys. Rev. B* **41**, 3905 (1990), https://doi.org/10.1103/PhysRevB.41.3905.
- Other nickel-centre table values: Nadolinny et al., *Crystals* **7**, 237 (2017), https://doi.org/10.3390/cryst7080237.

## Purpose

Constructs Spinach system and interaction specifications for the listed nickel-related diamond centres. This function selects a reported centre model and orientation; it does not run a spectrum or spin-dynamics simulation.

## Call and inputs

Call `[sys,inter]=diamond_ni(parameters)` with exactly one structure argument. Required fields are `parameters.centre` and `parameters.orientation`. Centre strings are case-insensitive and accepted values are `'w8'`, `'ne1'`–`'ne5'`, `'ne8'`, `'ab1'`–`'ab5'`, `'nol1'`, and `'nirim5'`. Orientation must be exactly `'111'`, `'110'`, or `'100'`; its crystal-plane normal is aligned with the applied field (`z`).

For W8, also supply `parameters.nickel` (a character isotope label or `'none'`) and `parameters.n_13c` (integer scalar 0–4). `'61Ni'` adds a `61Ni` nucleus with isotropic A = 0.65 mT; `'none'` adds no Ni nucleus. Another isotope string is placed in `sys.isotopes` without a hyperfine tensor; the function does not validate that label, so it must be a valid Spinach isotope name. For non-W8 centres, `n_13c` may be omitted or set to zero; a nonzero value is rejected.

Example: `[sys,inter]=diamond_ni(struct('centre','w8','orientation','111','nickel','none','n_13c',1));`

## Centre models

The code uses the following principal values. Hyperfine values listed in mT are converted to Hz internally using `abs(spin('E'))/(2*pi)*1e-3`. The g tensors are dimensionless.

| Centre | Electron / g principal values | Nuclear content |
|---|---|---|
| W8 | E4; isotropic g = 2.032; zero ZFS in the model | Optional 0–4 nearest-neighbour `13C`, A = [0.340, 0.340, 1.339] mT each; selected `61Ni` has A = 0.65 mT |
| NE1 | E; g = [2.1282, 2.0070, 2.0908] | two `14N` tensors, each [2.09, 1.43, 1.45] mT |
| NE2 | E; g = [2.1301, 2.0100, 2.0931] | three `14N` tensors: [2.10, 1.42, 1.41], [1.87, 1.18, 1.25], [0.18, 0.35, 0.25] mT |
| NE3 | E; g = [2.0729, 2.0100, 2.0476] | three `14N` tensors: [1.60, 1.24, 1.15], [0.66, 0.50, 0.50] twice, mT |
| NE4 | E; g = [2.0988, 2.0988, 2.0227] | no nuclear tensors are added by this case |
| NE5 | E; g = [2.0329, 2.0898, 2.0476] | two `14N` tensors, each [1.22, 0.98, 0.89] mT |
| NE8 | E; g = [2.0439, 2.1722, 2.0476] | four `14N` tensors, each [1.14, 0.78, 0.75] mT |
| AB1 | E; g = [2.0920, 2.0920, 2.0024] | no nuclear tensors are added by this case |
| AB2 | E; g = [2.0672, 2.0672, 2.0072] | no nuclear tensors are added by this case |
| AB3 | E; g = [2.1105, 2.0663, 2.0181] | no nuclear tensors are added by this case |
| AB4 | E; g = [2.0220, 2.0094, 2.0084] | no nuclear tensors are added by this case |
| AB5 | E3; g = [2.022, 2.022, 2.037]; axial ZFS parameter 1.132 T | no nuclear tensors are added by this case |
| NOL1 / NIRIM5 | E3; g = [2.002, 2.002, 2.0235]; axial ZFS parameter −6.10 T | no nuclear tensors are added by this case |

For the NE1/2/3/5/8 tensors, the principal-axis frame uses `alpha=14°` for NE1/2/3 and `27.5°` for NE5/8; AB3/AB4 use axes constructed from `[100]`, `[011]`, and `[0 −1 1]`; AB1/AB2, AB5, and NOL1/NIRIM5 use the source `frame_111`, while the NE centres use their alpha-defined frame. The review’s Table 2 labels its zero-field-splitting column in tesla, as the source comment notes. The NOL1/NIRIM5 and AB5 code multiplies the stated field-equivalent ZFS parameters by `abs(spin('E'))/(2*pi)` before constructing the tensor. Source comments connect these models to `D = −171 GHz` for NOL1/NIRIM5 (Nadolinny et al., *Diam. Relat. Mater.* **11**, 627 (2002)) and `D = 31.72 GHz` for AB5 (Landolt–Börnstein III/41A2a). The same comments note these splittings greatly exceed the 9.5 GHz X-band quantum; the AB5 splitting is below the 94 GHz W-band quantum of the shipped example.

## W8 carbon ordering and limits

The W8 carbons are placed, in order, on the `[111]`, `[1-1-1]`, `[-11-1]`, and `[-1-11]` nearest-neighbour bonds. Thus `n_13c` below four selects a particular isotopomer, not an average: at orientation `'111'`, the `[111]` carbon splits by 1.339 mT and each other listed carbon by 0.451 mT. To represent an isotopomer mixture, request all four and select or weight them in the calling script.

W8 is represented as a quartet with zero ZFS: the cited data in the source did not report W8 ZFS, and off-central transitions are treated as unresolved rather than explicitly modelled. The routine returns only its selected g, ZFS (where coded), and coupling matrices; it does not provide relaxation, line broadening, or mixture averaging.
