# interfaces/g2spinach.m

[Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/g2spinach.m) · [Spinach Wiki: g2spinach.m](https://spindynamics.org/wiki/index.php?title=g2spinach.m)

- Signature: `[sys,inter]=g2spinach(props,elements,references,options)`

## Inputs and selection

`props` is the property structure from `gparse` or `oparse`. `elements` is a cell array of element symbols, for example `{'H','N'}`. The routine selects matching atoms from `props.symbols` in geometry order. It takes their isotopes from `props.isotopes`, a cell array of isotope strings. Include `'E'` to select EPR import. The electron is appended as the final spin, not mapped to a geometry atom.

Supply `references` as a numeric vector with one value per entry of `elements`; these are absolute shielding values in ppm chosen to put each reference substance at zero chemical shift. Use reference calculations with the same method as the molecule. Zero values correspond to bare-nucleus referencing. `options` is optional; supported fields are `no_xyz`, `min_j`, `min_hfc`, and `purge` (the latter is relevant to EPR).

By default, selected standard-geometry rows are returned as per-spin coordinate cells in `inter.coordinates`; `options.no_xyz` suppresses them. The EPR electron coordinate is empty. This mapping assumes parser atom order and symbols are the desired atom identities; the routine does not use atom labels to disambiguate repeated elements.

## Ambiguous input

When the parser output is unambiguous the import proceeds; when it cannot be interpreted without a guess, the routine stops with an error that names the field of `props` to edit:

- `props.isotopes` is absent, is not a per-atom cell array, or is empty for a selected atom (CASTEP output from `c2spinach`, ORCA elements without a known default, elements above Z=103 in Gaussian): set the isotope strings.
- A selected atom has a non-magnetic isotope (Gaussian has no magnetic selenium and reports `80Se`): set `props.isotopes{n}` to the magnetic isotope to be imported.
- In NMR mode, `props.nqi{n}` is nonzero while the isotope of atom `n` has no quadrupole moment (ORCA labels hydrogens `1H` and prints deuterium quadrupolar tensors for them): set `props.isotopes{n}` to the isotope of the tensor, or empty `props.nqi{n}`.
- In NMR mode, printed J-couplings differ by more than 0.3% from those rebuilt from K-couplings with the stated isotopes, for pairs where both are nonzero (a checkpoint-read Gaussian job that carried isotopes the log does not state, an edited log): set `props.isotopes` to the isotopes of the calculation, or remove `props.j_couplings` to rebuild J-couplings for the stated isotopes. The scatter between the magnetic moment tables of Gaussian and Spinach is up to 0.14%.

Data that a mode ignores by design produce warnings: zero-field splitting is never imported; NMR mode ignores hyperfine couplings and the g-tensor; EPR mode ignores quadrupolar, scalar, and spin-rotation couplings.

## NMR import

When no electron is requested, an available shielding tensor `props.cst` becomes `inter.zeeman.matrix{n}=-props.cst{atom}+references(ref_index(n))I` in ppm, with the reference selected by element. For a nucleus with spin `I>1/2`, the reported quadrupolar tensor `props.nqi{atom}` in Hz is divided by `2I(2I-1)` and stored on the diagonal of `inter.coupling.matrix`. The converter does not use an electric field gradient to calculate a coupling for another isotope. Available spin-rotation tensors are copied to `inter.spinrot.matrix` in Hz as printed, for the isotopes of the calculation.

For scalar couplings, `k_couplings` takes precedence when present: the isotope-independent K values are converted to isotope-specific J values using both selected nuclei's gyromagnetic ratios and the source prefactor, and the code's factor of one-half. Otherwise `j_couplings` is selected directly with a factor of one-half. `options.min_j` retains only values whose absolute magnitude is strictly greater than the threshold. The returned scalar matrix is a cell array in Hz.

## EPR import

The electron's `props.g_tensor.matrix` is the only nonzero Zeeman tensor. Nuclear shifts and offsets are ignored. Each selected nucleus is coupled to the electron using its full hyperfine tensor for the isotope in the log. Half the tensor is converted from Gauss to Hz with `gauss2mhz` and stored in each direction. There is no isotope rescaling on import.

`options.min_hfc` clears each stored coupling whose Frobenius norm is below the threshold. With `options.purge='on'`, nuclei with no remaining coupling to the electron are removed from the isotope, Zeeman, coupling, and coordinate data. Without purge, these entries remain in the returned system.

## Outputs and reference values

`sys.isotopes` lists the selected isotope strings (and the appended electron in EPR mode). `inter` uses cells: `zeeman.matrix` holds one 3-by-3 tensor per spin; `coupling.matrix` is an N-by-N cell array of 3-by-3 tensors; `coupling.scalar` is an N-by-N cell array of scalar values; and coordinates and spin-rotation tensors are stored per spin. These fields are conditional on the applicable import path and available `props` fields. Scalar couplings and converted EPR hyperfine tensors are in Hz.

Inspect `sys.isotopes` before substituting an isotope. After import, `isoswap(sys,inter,spins,new_iso)` scales each selected spin's bilinear matrix, eigensystem, and scalar couplings by the new-to-old magnetogyric ratio. It rescales a nuclear quadrupolar coupling by the ratio of `Q/(2I(2I-1))` between the two isotopes, wipes with a warning the quadratic couplings it cannot transfer and all higher-order couplings, and refuses a source isotope with zero magnetogyric ratio.

For reference, the source documents these vacuum TMS absolute shielding values (13C, 1H):

| Method | 13C | 1H |
|---|---:|---:|
| GIAO, B3LYP/6-31G* | 189.6621 | 32.1833 |
| GIAO, B3LYP/6-311+G(2d,p) | 182.4485 | 31.8201 |
| GIAO, HF/6-31G* | 199.9711 | 32.5957 |
| GIAO, HF/6-311+G(2d,p) | 192.5828 | 32.0710 |
| CSGT, B3LYP/6-31G* | 188.5603 | 29.1952 |
| CSGT, B3LYP/6-311+G(2d,p) | 182.1386 | 31.7788 |
| CSGT, HF/6-31G* | 196.8670 | 29.5517 |
| CSGT, HF/6-311+G(2d,p) | 192.5701 | 31.5989 |

In EPR mode the reference values are ignored. The function returns `sys` and `inter`; it does not run an electronic-structure calculation. See the linked MATLAB source for the exact input guards and transformations.
