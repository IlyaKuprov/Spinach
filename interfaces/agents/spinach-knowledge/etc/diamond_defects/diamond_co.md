# etc/diamond_defects/diamond_co.m

- Signature: `[sys,inter]=diamond_co(parameters)`
- [Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=diamond_co.m)
- Magnetic parameters: [Nadolinny et al., *Crystals* 7, 237 (2017)](https://doi.org/10.3390/cryst7080237).

## Purpose

Build a Spinach spin-system and interaction specification for a cobalt-related defect in diamond, using the tabulated magnetic parameters for either the `o4` or `nlo2` cobalt centre.

## Input

`parameters` is a structure with two required character fields:

- `centre` — `o4` or `nlo2` (case-insensitive).
- `orientation` — `111`, `110`, or `100`, selecting the crystal-plane normal aligned with the magnetic field.

## Construction

The routine selects the centre-specific principal values for the electron g tensor and the 59Co hyperfine tensor, constructs the cobalt principal-axis frame, and rotates it to the requested field orientation.

## Outputs

- `sys` — Spinach system specification for an electron and a 59Co nucleus.
- `inter` — Spinach interaction specification containing their Zeeman and hyperfine coupling tensors.
