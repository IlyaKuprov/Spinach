# etc/diamond_defects/diamond_gev0.m

- MATLAB implementation: [etc/diamond_defects/diamond_gev0.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/diamond_defects/diamond_gev0.m)

- Signature: `[sys,inter]=diamond_gev0(parameters)`
- [Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=diamond_gev0.m)
- Magnetic parameters: Nadolinny et al., *Phys. Status Solidi A* **213**, 2623 (2016), [https://doi.org/10.1002/pssa.201600211](https://doi.org/10.1002/pssa.201600211)
- Source author: alexey.bogdanov@weizmann.ac.il

## Purpose

Build Spinach system and interaction structures for the GeV0 defect in diamond, including the electron Zeeman and zero-field-splitting terms and, when selected, a germanium nucleus.

## Call and inputs

`[sys,inter]=diamond_gev0(parameters)` takes exactly one structure containing:

- `parameters.germanium` — character string. Use `'73Ge'` for the parameterised isotope, `'none'` to omit a nucleus, or another germanium isotope label. The code passes any other non-`'none'` string through as the isotope label; it only assigns a hyperfine tensor for exact `'73Ge'`.
- `parameters.orientation` — exact character string `'111'`, `'110'`, or `'100'`; the chosen crystal plane normal is aligned with the magnetic field.

## Magnetic tensors and units

The electron isotope label is `'E3'`. Its principal `g` values are [2.0027, 2.0027, 2.0025]. The source constructs the zero-field-splitting matrix with `zfs2mat(80.3·hz_per_mt,0,0,0,0)`, where `hz_per_mt=abs(spin('E'))/(2π)·10⁻³`; the 80.3 value is thus converted from mT to Spinach frequency units before assembly. The fixed electron frame used to build the tensors is `[-1/√2, -1/√6, 1/√3; 1/√2, -1/√6, 1/√3; 0, 2/√6, 1/√3]`; it is then transformed into the selected crystallographic orientation. The rotated ZFS matrix is converted with `mat2ias` for the interaction structure.

For exact `'73Ge'`, an isotropic hyperfine tensor `1.64·hz_per_mt·I` is added (1.64 mT expressed in Spinach frequency units). For another isotope label no hyperfine tensor is set by this source.

## Outputs

- `sys` — electron isotope `'E3'`, followed by the optional germanium isotope unless `germanium='none'`.
- `inter` — rotated electron Zeeman tensor at `inter.zeeman.matrix{1}`; ZFS interaction at `inter.coupling.matrix{1,1}`; and, only for `'73Ge'`, electron–germanium hyperfine tensor at `inter.coupling.matrix{1,2}`.

The routine checks that both fields are character strings and that orientation is one of the three listed values. It does not verify that a non-special germanium string names a supported isotope.
