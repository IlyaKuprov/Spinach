# etc/diamond_defects/diamond_n2vm.m

- Signature: `[sys,inter]=diamond_n2vm(parameters)`
- [Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=diamond_n2vm.m)
- Magnetic parameters: [Green et al., *Physical Review B* 92, 165204 (2015)](https://doi.org/10.1103/PhysRevB.92.165204).

## Purpose

Construct a Spinach model of the N2V− defect in diamond, with two nitrogen nuclei and optional reported nearest-neighbour 13C couplings.

## Inputs

`parameters` is a structure with three required fields:

- `nitrogen` — `14N` or `15N`, used for both nitrogen nuclei.
- `orientation` — `111`, `110`, or `100`, selecting the crystal-plane normal aligned with the magnetic field.
- `include_13c` — logical scalar; `true` adds the two reported 13C nuclei and their hyperfine couplings.

## Construction

The routine assembles the electron g tensor and orientation-dependent nitrogen hyperfine tensors. For 14N it also includes the nuclear quadrupole tensors. The selected crystal orientation is aligned to the field before the interactions are placed in the Spinach structures.

## Outputs

- `sys` — Spinach system specification for the electron, two nitrogen nuclei, and any included 13C nuclei.
- `inter` — Spinach interaction specification with the electron Zeeman tensor and the applicable hyperfine and quadrupole couplings.
