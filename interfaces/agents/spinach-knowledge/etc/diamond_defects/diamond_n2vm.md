# etc/diamond_defects/diamond_n2vm.m

- MATLAB implementation: [etc/diamond_defects/diamond_n2vm.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/diamond_defects/diamond_n2vm.m)

- Signature: `[sys,inter]=diamond_n2vm(parameters)`
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=diamond_n2vm.m)
- Magnetic parameters: Green et al., *Phys. Rev. B* **92**, 165204 (2015), https://doi.org/10.1103/PhysRevB.92.165204.

## Purpose

Constructs Spinach `sys` and `inter` specifications for the N2V`−` diamond defect: one electron, two nitrogen nuclei, and optionally the two reported nearest-neighbour `13C` couplings. It returns model specifications, not a Hamiltonian simulation or a pulse sequence.

## Call and inputs

Call `[sys,inter]=diamond_n2vm(parameters)` with exactly one structure argument. Required fields:

- `parameters.nitrogen`: exactly `'14N'` or `'15N'`.
- `parameters.orientation`: `'111'`, `'110'`, or `'100'`; the corresponding crystal-plane normal is aligned with the applied field (laboratory `z` axis).
- `parameters.include_13c`: logical scalar `true` or `false`.

Example: `[sys,inter]=diamond_n2vm(struct('nitrogen','14N','orientation','111','include_13c',false));`

## Model and interactions

The electron g principal values are `[2.00345, 2.00274, 2.00271]`. Each nitrogen has a hyperfine tensor with principal values `[3.47, 4.51, 4.09] MHz`; the second is rotated by 180° about the crystal `z` axis. For `14N` the code scales these hyperfine values by `spin('14N')/spin('15N')` and adds the reported quadrupole tensor with principal value `−5.0 MHz` in its specified frame. For `15N` no quadrupole tensor is added. The nitrogen frame includes a −3.5° rotation about `[1 −1 0]`; the two sites differ by a 180° rotation about the crystal z axis.

When `include_13c` is true, two `13C` nuclei are appended. Their hyperfine principal values are `[202.3, 202.3, 317.5] MHz`, with the second tensor rotated by 180° about `z`. Its principal-axis frame includes a +2.0° rotation about `[-1 −1 0]`. The routine constructs the crystal/principal-axis frames, rotates the electron Zeeman and nuclear couplings into the selected field orientation, and stores them in Spinach's interaction matrices.

## Outputs and scope

- `sys` contains the electron and selected nuclear isotope labels.
- `inter` contains the electron Zeeman tensor plus electron–nuclear hyperfine and, for `14N`, quadrupole information.

The number and isotope of nuclei are fixed by these parameters; this is not a general multi-carbon builder.
