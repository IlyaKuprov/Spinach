# etc/diamond_defects/diamond_gev0.m

- Signature: `[sys,inter]=diamond_gev0(parameters)`
- [Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=diamond_gev0.m)
- Magnetic parameters: [Nadolinny et al., *Physica Status Solidi A* 213, 2623 (2016)](https://doi.org/10.1002/pssa.201600211).

## Purpose

Construct a Spinach model of the GeV0 defect in diamond, including the electron g tensor and zero-field splitting, with an optional germanium nuclear spin.

## Input

`parameters` is a structure with two required character fields:

- `germanium` — `73Ge`, `none`, or another germanium isotope name. `73Ge` receives the tabulated hyperfine coupling; another non-`none` value adds a nucleus with that isotope label.
- `orientation` — `111`, `110`, or `100`, selecting the crystal-plane normal aligned with the magnetic field.

## Construction

The routine rotates the electron g and zero-field-splitting tensors to the selected orientation. When a germanium nucleus is included, it also sets its hyperfine coupling.

## Outputs

- `sys` — Spinach system specification containing the electron and, when requested, the germanium nucleus.
- `inter` — Spinach interaction specification containing the electron Zeeman tensor, zero-field splitting, and any applicable hyperfine coupling.
