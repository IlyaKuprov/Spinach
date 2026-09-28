# examples/dnp_sol/solid_effect_field_scan_1.m

- Signature: `solid_effect_field_scan_1()`

## Purpose

Sweeps magnetic-field offset for a gadolinium-containing DNP system and computes the steady-state (^{15}mathrm{N}) (S_z) expectation value. The source estimates the calculation time as minutes.

## Physical / mathematical content

- The spin system comprises the Spinach `E8` gadolinium electron spin and a (^{15}mathrm{N}) nucleus, with nominal field 9.4509 T and temperature 40.2 K.
- The electron Zeeman principal values are all 1.9918. The source specifies a 570 MHz axial electron–nuclear coupling tensor and places the spins 3.00 Å apart along (x).
- Relaxation uses `t1_t2`, diagonal terms retained, zero equilibrium, and the listed rates `r1_rates={1e4 1e1}` and `r2_rates={1e7 1e3}`.
- The 512-point scan spans field offsets from -0.08 to +0.08 T. The detected signal is the nitrogen (S_z) expectation value.

## Numerical / algorithmic content

Builds the spin system with a full `sphten-liouv` basis and the listed `E8` projections, then calls `powder` with the `dnp_field_scan` callback in ESR mode. The powder grid is `rep_2ang_6400pts_sph`; the requested method is `backslash`, and the calculation requests an anisotropic-equilibrium term.

## Implementation structure

The function sets the field, isotopes, electron Zeeman and coupling tensors, coordinates, basis and relaxation model; creates the Spinach system; configures the electron microwave operators, frequency, power, field-offset points, powder grid and nitrogen detection state; runs the powder field scan; and plots the real nitrogen signal against field offset.
