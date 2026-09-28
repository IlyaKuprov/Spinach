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

The routine selects the centre-specific principal values for the electron g tensor and the 59Co hyperfine tensor, constructs the cobalt principal-axis frame, and rotates it to the requested field orientation. The values encoded in the routine are:

| Centre | Principal g values | 59Co hyperfine principal values (mT) | Frame angle α (degrees) |
| --- | --- | --- | --- |
| `o4` | 2.3463, 1.8438, 1.7045 | 8.86, 6.43, 5.82 | 29 |
| `nlo2` | 2.3277, 1.7982, 1.7149 | 8.24, 6.57, 5.76 | 28 |

The hyperfine values are converted from millitesla to the frequency units used by Spinach; `α` sets the cobalt principal-axis frame.

## Outputs

- `sys` — Spinach system specification for an electron and a 59Co nucleus.
- `inter` — Spinach interaction specification containing their Zeeman and hyperfine coupling tensors.
