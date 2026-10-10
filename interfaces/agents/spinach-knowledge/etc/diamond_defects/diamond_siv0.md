# etc/diamond_defects/diamond_siv0.m

- MATLAB implementation: [etc/diamond_defects/diamond_siv0.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/diamond_defects/diamond_siv0.m)

**Call:** `[sys,inter] = diamond_siv0(parameters)`

Builds the SiV0 split-vacancy centre model for diamond: one E3 electronic-spin site, optional silicon, and a selectable number of nearest-neighbour `13C` spins.

## Inputs

- `parameters.silicon`: `'29Si'` adds the explicitly parameterised silicon spin; `'none'` omits it; another character isotope label adds that isotope without a silicon hyperfine tensor from this routine.
- `parameters.orientation`: `'111'`, `'110'` or `'100'`; the chosen crystal direction is aligned with the magnetic field.
- `parameters.n_13c`: integer from 0 through 6, selecting the number of reported nearest-neighbour carbon couplings.

All three fields are required. The source validates the orientation against those exact character values.

## Model details

The electron label is `E3`; the principal `g` values are `[2.0035, 2.0035, 2.0042]`, and the axial ZFS is passed to `zfs2mat` as `1000e6`. For `'29Si'`, the source adds an anisotropic tensor with principal values `[78.9e6, 78.9e6, 76.3e6]`. The source does not annotate units for these numerical arguments, so they are reported here as coded rather than relabelled.

For requested carbons, each tensor has coded principal values `[30.2e6, 30.2e6, 66.2e6]`, oriented along one of three split-vacancy dangling-bond directions. Selection cycles through those three directions with `mod(n-1,3)+1`; counts 4–6 therefore reuse the direction sequence. Orientation rotations are applied to the electron tensors and the resulting Zeeman/ZFS matrices are placed in `inter`.

## Outputs and limits

- `sys`: electron and whichever silicon/carbon isotopes were requested.
- `inter`: Zeeman, ZFS and available nuclear hyperfine matrices.

The routine's silicon special case is specifically `'29Si'`; other isotope labels do not inherit its explicit hyperfine values. It constructs interaction tensors, not a relaxed coordinate structure.

## Reference

Magnetic parameters are attributed to Edmonds et al., *Physical Review B* **77**, 245205 (2008), [doi:10.1103/PhysRevB.77.245205](https://doi.org/10.1103/PhysRevB.77.245205).

[Spin Dynamics Wiki source page](https://spindynamics.org/wiki/index.php?title=diamond_siv0.m).
