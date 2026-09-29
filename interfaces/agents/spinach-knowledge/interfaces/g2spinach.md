# interfaces/g2spinach.m

[Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/g2spinach.m) · [Spinach Wiki: g2spinach.m](https://spindynamics.org/wiki/index.php?title=g2spinach.m)

- Signature: `[sys,inter]=g2spinach(props,particles,references,options)`

## Inputs and selection

`props` is the property structure from `gparse` or `oparse`. `particles` is a cell array of two-entry cells, each containing an element symbol and a Spinach isotope string, for example `{{'H','1H'},{'N','15N'}}`. The routine walks `props.symbols` in atom order and selects entries whose element symbol matches a requested element; it returns the corresponding requested isotope in `sys.isotopes`. Add `{'E','E'}` to select EPR import. In EPR mode the electron is appended as the final spin, not mapped to a geometry atom.

Supply `references` as a numeric vector with one value per entry of `particles`; these are absolute shielding values in ppm chosen to put each reference substance at zero chemical shift. Use reference calculations with the same method as the molecule. Zero values correspond to bare-nucleus referencing. `options` is optional; supported fields are `no_xyz`, `min_j`, `min_hfc`, and `purge` (the latter is relevant to EPR). The implementation's input checks require `props` to be a structure, `particles` to be a cell array, and the numeric reference vector to have matching length.

By default, selected standard-geometry rows are returned as per-spin coordinate cells in `inter.coordinates`; `options.no_xyz` suppresses them. The EPR electron coordinate is empty. This mapping assumes parser atom order and symbols are the desired atom identities; the routine does not use atom labels to disambiguate repeated elements.

## NMR import

When no electron is requested, an available shielding tensor `props.cst` becomes `inter.zeeman.matrix{n}=-props.cst{atom}+references(ref_index(n))I` in ppm, with the matching reference selected from the particle list. Available nuclear quadrupole tensors are placed on the diagonal of `inter.coupling.matrix`; for nuclei with spin `I>1/2`, the parser tensor is divided by `2I(2I-1)`. Available spin-rotation tensors are copied to `inter.spinrot.matrix`. These parser fields use Hz for quadrupolar and spin-rotation tensors.

For scalar couplings, `k_couplings` takes precedence when present: the isotope-independent K values are converted to isotope-specific J values using both selected nuclei's gyromagnetic ratios and the source prefactor, and the code's factor of one-half. Otherwise `j_couplings` is selected directly with a factor of one-half; this branch prints a warning that those values are isotope-specific and are not rescaled. `options.min_j` retains only values whose absolute magnitude is strictly greater than the threshold. The returned scalar matrix is a cell array in Hz.

## EPR import

The electron's `props.g_tensor.matrix` is the only nonzero Zeeman tensor. Nuclear shifts and offsets are ignored. Each selected nucleus is coupled symmetrically to the electron using its full hyperfine tensor: the source Gauss tensor is converted to Hz with `gauss2mhz`, then scaled by target/source gyromagnetic ratio. Nonempty HFCs require an explicit source isotope in `props.isotopes` (Gaussian mass numbers or ORCA isotope strings); missing, invalid, or zero-gyromagnetic-ratio sources are rejected. The converter uses Spinach's `spin` and `gauss2mhz` routines for isotope and field-unit conversions.

`options.min_hfc` clears each coupling whose scaled tensor's Frobenius norm is below the threshold. With `options.purge='on'`, nuclei with no remaining coupling to the electron are removed from the isotope, Zeeman, coupling, and coordinate data. Without purge, these entries remain in the returned system.

## Outputs and reference values

`sys.isotopes` lists the selected isotope strings (and the appended electron in EPR mode). `inter` uses cells: `zeeman.matrix` holds one 3-by-3 tensor per spin; `coupling.matrix` is an N-by-N cell array of 3-by-3 tensors; `coupling.scalar` is an N-by-N cell array of scalar values; and coordinates and spin-rotation tensors are stored per spin. These fields are conditional on the applicable import path and available `props` fields. Scalar couplings and converted EPR hyperfine tensors are in Hz.

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
