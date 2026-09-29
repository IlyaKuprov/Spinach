# etc/diamond_defects/diamond_p.m

- MATLAB implementation: [etc/diamond_defects/diamond_p.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/diamond_defects/diamond_p.m)

## Purpose and call

`[sys,inter]=diamond_p(parameters)` constructs Spinach spin-system and interaction specifications for the phosphorus-related diamond centres described by Nadolinny et al., *Crystals* **7**, 237 (2017), <https://doi.org/10.3390/cryst7080237>.

## Inputs and constraints

`parameters` must be one structure with these fields:

- `centre`: character string naming `'ma1'`, `'np1'`, `'np2'`, `'np3'`, `'np4'`, `'np5'`, `'np6'`, `'np8'`, or `'np9'`; matching is case-insensitive.
- `orientation`: character string `'111'`, `'110'`, or `'100'`, specifying the crystal-plane normal aligned with the magnetic field.
- `include_13c`: scalar logical controlling the reported MA1 13C hyperfine coupling. It defaults to `false` if absent; `true` is accepted only for MA1.

## Model and units

The centre selects an electron g tensor and nuclear hyperfine tensors. Each centre includes 31P; NP1–NP3 also include 14N; NP8 includes two 31P nuclei; MA1 can include 13C when requested. Tabulated hyperfine principal values are converted from mT to frequency units using `abs(spin('E'))/(2*pi)*1e-3`. Tensors are placed in the frames specified for each centre and rotated to align the selected crystal-plane normal with the field. The routine populates the electron Zeeman matrix and electron–nuclear coupling matrices.

## Outputs and limitations

- `sys`: Spinach system specification structure containing the electron and selected nuclear isotopes.
- `inter`: Spinach interaction specification structure containing the electron Zeeman tensor and electron–nuclear hyperfine couplings.

No zero-field-splitting or nuclear quadrupole interaction is assigned. The function returns specifications, not a simulated spectrum. Source documentation: <https://spindynamics.org/wiki/index.php?title=diamond_p.m>.

**Frame clarification:** NP6 uses the trigonal frame for its g tensor but the identity frame for its 31P hyperfine tensor; do not assume all tensors for one centre share a frame.
