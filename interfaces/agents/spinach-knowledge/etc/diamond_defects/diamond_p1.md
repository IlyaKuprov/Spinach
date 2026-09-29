# etc/diamond_defects/diamond_p1.m

- MATLAB implementation: [etc/diamond_defects/diamond_p1.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/diamond_defects/diamond_p1.m)

## Purpose and call

Constructs a P1-centre spin system for diamond:

`[sys,inter]=diamond_p1(parameters)`

`parameters` must be a structure. Its supported fields are:

- `orientation`: `'111'`, `'110'`, or `'100'`; the specified crystal-plane normal is aligned with the magnetic field. Defaults to `'111'`.
- `nitrogen`: `'14N'` or `'15N'`. Defaults to `'14N'`.

## Physical and numerical model

The system contains an electron (`'E'`) and the selected nitrogen isotope. The electron g-tensor principal values in the trigonal frame are 2.00220, 2.00220, and 2.00218. For `'14N'`, the electron–nitrogen hyperfine values are 81.3, 81.3, and 114.0 MHz, and a nitrogen quadrupolar tensor is assigned as `zfs2mat(-3.97e6,0,0,0,0)`. For `'15N'`, the hyperfine values are −114.0, −114.0, and −159.9 MHz; no nitrogen quadrupolar tensor is assigned. Tensors are transformed from the trigonal principal-axis frame with a rotation aligning the selected crystal direction with `[0 0 1]`.

Magnetic parameters: Nir-Arad et al., *Phys. Chem. Chem. Phys.* **26**, 27633 (2024), <https://doi.org/10.1039/d4cp03055a>; Smith et al., *Phys. Rev.* **115**, 1546 (1959), <https://doi.org/10.1103/PhysRev.115.1546>.

## Outputs and limitations

- `sys`: Spinach system specification structure.
- `inter`: Spinach interaction specification structure containing the electron Zeeman tensor and isotope-dependent coupling tensors.

The function builds specifications, not a simulation, and does not set a magnetic-field strength. Although both fields default when omitted from the structure, the `parameters` structure argument itself is still required. Supplied `orientation` must be a character array; supplied `nitrogen` is not separately type-checked. Unsupported orientation and isotope values error. Source documentation: <https://spindynamics.org/wiki/index.php?title=diamond_p1.m>.

**Orientation clarification:** `'111'`, `'110'`, and `'100'` identify the crystal-plane normal aligned with the field; they are not field-strength settings.
