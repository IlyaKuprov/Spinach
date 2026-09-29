# etc/diamond_defects/diamond_n_inter.m

- MATLAB implementation: [etc/diamond_defects/diamond_n_inter.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/diamond_defects/diamond_n_inter.m)

- Signature: `[sys,inter]=diamond_n_inter(parameters)`
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=diamond_n_inter.m)
- Magnetic parameters: Felton et al., *J. Phys.: Condens. Matter* **21**, 364212 (2009), https://doi.org/10.1088/0953-8984/21/36/364212.

## Purpose

Builds the electron–nitrogen spin-system specifications for the War9 or War10 nitrogen-interstitial centre in diamond. It returns `sys` and `inter`; it does not perform a spin-dynamics calculation.

## Call and inputs

Call `[sys,inter]=diamond_n_inter(parameters)` with exactly one structure argument. Required fields:

- `parameters.centre`: character string `'war9'` or `'war10'` (case-insensitive; the function lowercases this field).
- `parameters.orientation`: exactly `'111'`, `'110'`, or `'100'`, specifying the crystal-plane normal aligned with the applied field (`z`).
- `parameters.nitrogen`: exactly `'14N'` or `'15N'`.

Example: `[sys,inter]=diamond_n_inter(struct('centre','war9','orientation','111','nitrogen','14N'));`

## Centre tensors and isotope treatment

The tabulated g and hyperfine tensors are expressed in centre frames and rotated into the requested field orientation. Principal values encoded by the source are:

| Centre | g principal values | Nitrogen A principal values |
|---|---|---|
| War9 | `[2.00343, 2.00272, 2.00268]` | `[8.30, 7.85, 8.17] MHz` |
| War10 | `[2.00344, 2.00272, 2.00269]` | `[1.00, −1.01, 0.00] MHz` |

For `15N` these A values are used directly. For `14N` the hyperfine tensor is multiplied by `spin('14N')/spin('15N')`. War9 uses a shared g/hyperfine frame built from directions `(θ,φ)=(90°,45°),(180°,45°),(90°,315°)`. For War10 the g tensor stays in that frame, while the hyperfine frame uses `(44.8°,45.0°)`, `(134.8°,45.0°)`, and `(90°,315°)`. Although `14N` is quadrupolar, this routine does not add a nuclear quadrupole interaction: its returned coupling matrix contains only the electron–nitrogen hyperfine tensor.

## Outputs and scope

- `sys.isotopes` identifies the electron (`'E'`) and selected nitrogen isotope.
- `inter.zeeman.matrix{1}` holds the rotated electron g tensor; `inter.coupling.matrix{1,2}` holds the rotated A tensor.

The War10 tensor is anisotropic and has a negative second principal value; it should not be read as an isotropic positive coupling. The code creates no additional nuclei or NQI terms.
