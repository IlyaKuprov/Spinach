# etc/diamond_defects/diamond_r2.m

- MATLAB implementation: [etc/diamond_defects/diamond_r2.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/diamond_defects/diamond_r2.m)

**Call:** `[sys,inter] = diamond_r2(parameters)`

Builds the R2 self-interstitial model as one E3 electronic-spin site in Spinach.

## Inputs

- `parameters.orientation`: `'111'`, `'110'` or `'100'`; the selected crystal direction is aligned with the magnetic field.
- `parameters.d_sign`: real numeric scalar multiplying the axial zero-field-splitting parameter. The source validates scalar reality, but does not restrict this multiplier to exactly `+1` or `-1`.

Both fields are required. Orientation values are matched exactly.

## Model and outputs

The source assigns the electron isotope label `E3`, constructs principal g values `[2.0019, 2.0019, 2.0021]`, and forms an axial zero-field-splitting tensor with `zfs2mat(parameters.d_sign * 4173e6, 0, 0, 0, 0)`. The numeric `4173e6` argument is passed directly to `zfs2mat`; the routine's comments do not independently state its unit. A rotation maps the chosen crystal direction onto the field axis, then rotates both electron tensors into that frame.

- `sys`: the single-site electronic-spin specification (`E3`).
- `inter`: electron Zeeman tensor and the axial ZFS interaction, with ZFS converted through `mat2ias` for Spinach's interaction structure.

This helper does not add nuclei or expose a magnitude range for `d_sign`; callers can choose any real scalar accepted by the input check.

## Reference

Magnetic parameters are attributed in the source to Hunt et al., *Physical Review B* **61**, 3863 (2000), [doi:10.1103/PhysRevB.61.3863](https://doi.org/10.1103/PhysRevB.61.3863).

[Spin Dynamics Wiki source page](https://spindynamics.org/wiki/index.php?title=diamond_r2.m).
