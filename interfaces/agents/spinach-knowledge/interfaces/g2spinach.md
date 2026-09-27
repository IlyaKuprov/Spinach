# interfaces/g2spinach.m

- Signature: `[sys,inter]=g2spinach(props,particles,references,options)`

## Purpose

Converts parsed Gaussian or ORCA electronic-structure properties into Spinach isotope and interaction data. Including an electron in `particles` selects EPR import; otherwise the routine imports NMR parameters.

## Import behavior

For NMR, nuclear chemical-shift tensors are formed from the negative shielding tensor plus the supplied reference value. Available quadrupolar tensors and spin-rotation tensors are imported. Isotope-independent `k_couplings` are converted to isotope-specific scalar couplings using the spins' gyromagnetic ratios; supplied `j_couplings` are used as isotope-specific values and trigger a warning.

For EPR, the electron is appended to the isotope list, its g-tensor is imported, and the selected nuclei's hyperfine tensors are converted to Hz and scaled by the target/source nuclear gyromagnetic-ratio. This scaling is applied before thresholding or purging. Each selected nucleus with a nonempty source hyperfine tensor must have explicit source-isotope metadata in `props.isotopes` (Gaussian mass numbers or ORCA isotope strings); malformed or missing metadata and zero-gamma source isotopes are rejected. The EPR branch sets nuclear chemical shifts to zero and ignores the reference values and offset specification. Molecular coordinates come from `props.std_geom`; `options.no_xyz=1` suppresses them, and the electron has no molecular coordinate.

## Parameters / inputs

- `props`: parsed output from `gparse` or `oparse`, with properties needed for the selected mode. EPR import requires explicit source-isotope metadata for every selected nonempty hyperfine tensor.
- `particles`: cell array of element/isotope pairs to import, for example `{{'H','1H'},{'N','15N'}}`. Including `{'E','E'}` selects EPR mode.
- `references`: vector of absolute shielding values for the reference substances, one per particle, to place at zero ppm; calculate these with the same electronic-structure method. References are ignored in EPR mode.

  The source gives the following absolute isotropic shielding values for tetramethylsilane in vacuum:

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

  This reference setting is ignored when electrons are present.
- `options.min_j`: scalar-coupling threshold in Hz; NMR couplings with absolute value at or below this threshold are zeroed.
- `options.min_hfc`: EPR hyperfine threshold in Hz; tensors whose Frobenius norm is below it are removed.
- `options.purge`: when set to `'on'` in EPR mode, removes nuclear spins with no remaining hyperfine coupling after thresholding.
- `options.no_xyz`: when set to 1, omits coordinate information while retaining interaction tensors.

## Outputs

- `sys.isotopes`: imported isotope labels.
- `inter.coordinates`: imported molecular coordinates in Angstrom when enabled; the EPR electron has an empty coordinate entry.
- `inter.zeeman.matrix`: one 3-by-3 matrix per spin, in ppm for nuclei and as the g-tensor for the electron.
- `inter.coupling.matrix`: 3-by-3 coupling tensors in Hz, including imported hyperfine or quadrupolar interactions.
- `inter.coupling.scalar`: scalar couplings in Hz when supplied.
- `inter.spinrot.matrix`: imported spin-rotation tensors when present.

## Source

[Spinach Wiki: g2spinach.m](https://spindynamics.org/wiki/index.php?title=g2spinach.m)
