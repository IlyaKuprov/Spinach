# etc/diamond_defects/diamond_nvm_gs.m

- MATLAB implementation: [etc/diamond_defects/diamond_nvm_gs.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/diamond_defects/diamond_nvm_gs.m)

## Purpose and signature

Constructs a diamond NV-centre ground-state spin system using magnetic parameters from S. Felton et al., *Phys. Rev. B* **79**, 075203 (2009), <https://doi.org/10.1103/PhysRevB.79.075203>.

`[sys,inter]=diamond_nvm_gs(parameters)`

## Parameters and constraints

`parameters` must be a structure; both fields are optional:

- `parameters.orientation`: `'111'`, `'110'`, or `'100'`, specifying the crystal-plane normal aligned with the magnetic field. Defaults to `'111'`. If supplied, it must be a character string; an unsupported value errors.
- `parameters.nitrogen`: `'14N'` or `'15N'`; defaults to `'14N'`. Unsupported values error; the source does not separately check this field's type.

## Physical and numerical model

The isotope list is `{'E3',parameters.nitrogen}`. Tensors are defined in a trigonal principal-axis frame and rotated so the selected crystal direction aligns with the field. The electron g-tensor principal values are 2.0031, 2.0031, and 2.0029; the electron zero-field-splitting tensor is `zfs2mat(2872e6,0,0,0,0)` (D = 2872 MHz). For `'14N'`, the hyperfine principal values are −2.70, −2.70, and −2.14 MHz, with a quadrupolar tensor `zfs2mat(-5.01e6,0,0,0,0)` (D = −5.01 MHz). For `'15N'`, the hyperfine values are 3.65, 3.65, and 3.03 MHz; no nitrogen quadrupolar tensor is assigned.

## Outputs and scope

- `sys`: Spinach system specification structure containing the isotope list.
- `inter`: Spinach interaction specification structure containing the electron Zeeman tensor and the zero-field-splitting and isotope-dependent coupling tensors.

The function builds specifications; it does not calculate a spectrum or set a magnetic-field magnitude. Only the listed orientations and isotopes have implemented branches. Source documentation: <https://spindynamics.org/wiki/index.php?title=diamond_nvm_gs.m>.
