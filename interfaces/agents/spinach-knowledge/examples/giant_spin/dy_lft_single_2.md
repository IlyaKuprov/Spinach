# examples/giant_spin/dy_lft_single_2.m

- MATLAB implementation: [examples/giant_spin/dy_lft_single_2.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/giant_spin/dy_lft_single_2.m)

- Source: [examples/giant_spin/dy_lft_single_2.m](../../../../../examples/giant_spin/dy_lft_single_2.m)
- Signature: `dy_lft_single_2()` (no input or output arguments)

## Purpose

Demonstrates the zero-field-splitting (ZFS) limit for relaxation theory in a single-lanthanide model; the source describes it as a figure from forthcoming papers. The source estimates calculation time in hours; it was not timed here.

## Spin model and construction

The system is one `E16` Dy ion with `sys.magnet=1.0`. A real anisotropic g tensor is built from `D=[1.322766699, 1.324261429, 1.328750739]` and the source rotation matrix `V`. The MOLCAS ligand-field arrays at ranks 2, 4, and 6 are converted with `icm2hz` and `stev2sph`, rotated using the Euler angles derived from `R`, and supplied through `inter.giant.coeff`; the associated giant-spin Euler angles are zero.

The basis uses `formalism='zeeman-hilb'` and `approximation='none'`, followed by `create` and `basis`.

## Field scan

The script calls `fieldscan_enlev` with `fields=[0 500]`, `npoints=1000`, `orientation=[0 0 0]`, and `nstates=16`. The called `fieldscan_enlev` function interprets the field endpoints in tesla and the orientation Euler angles in radians. This is an energy-level scan: the function has no output arguments, and the call does not capture or plot a returned spectrum. The MOLCAS ligand-field arrays are passed to `icm2hz` as inverse-centimetre values; that function converts them to Hz.
