# etc/diamond_defects/diamond_nv0_es.m

- MATLAB implementation: [etc/diamond_defects/diamond_nv0_es.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/diamond_defects/diamond_nv0_es.m)

- Signature: `[sys,inter]=diamond_nv0_es(parameters)`
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=diamond_nv0_es.m)
- Magnetic parameters: Felton et al., *Phys. Rev. B* **77**, 081201 (2008), https://doi.org/10.1103/PhysRevB.77.081201.

## Purpose

Builds the NV`0` excited-state electron–nitrogen spin-system specifications. The effective electron is a quartet (`E4`, spin 3/2), coupled to one nitrogen. The function returns model matrices; it does not propagate dynamics or simulate an experiment.

## Call and inputs

Call `[sys,inter]=diamond_nv0_es(parameters)` with exactly one structure argument. Required fields:

- `parameters.orientation`: exactly `'111'`, `'110'`, or `'100'`; the corresponding crystal-plane normal is aligned with the applied field (`z`).
- `parameters.nitrogen`: exactly `'14N'` or `'15N'`.

Example: `[sys,inter]=diamond_nv0_es(struct('orientation','111','nitrogen','14N'));`

## Magnetic model

In the principal frame the electron g values are `[2.0035, 2.0035, 2.0029]`, axial zero-field splitting is `D=1685 MHz`, and the `15N` hyperfine values are `[−23.8, −23.8, −35.7] MHz`. The `14N` hyperfine tensor is obtained by multiplying those values by `spin('14N')/spin('15N')`. No `14N` nuclear quadrupole interaction is added because the source says none is reported for this state.

The selected plane normal is rotated to the laboratory field axis. The output specifies electron and nitrogen isotope labels, the rotated electron Zeeman tensor, electron zero-field splitting converted to Spinach's interaction representation, and electron–nitrogen hyperfine coupling.

## Outputs and scope

- `sys` contains `'E4'` and the selected nitrogen isotope.
- `inter` contains the Zeeman and zero-field-splitting electron terms and the nitrogen hyperfine tensor.

This is specifically the NV`0` excited-state parameterisation; it is not the NV`−` ground-state model and includes no NQI term.
