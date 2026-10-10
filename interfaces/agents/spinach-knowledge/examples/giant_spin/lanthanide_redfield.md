# examples/giant_spin/lanthanide_redfield.m

- MATLAB implementation: [examples/giant_spin/lanthanide_redfield.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/giant_spin/lanthanide_redfield.m)

- Signature: `lanthanide_redfield()`
- Source: `examples/giant_spin/lanthanide_redfield.m`

## Model and parameters

This example treats Gd(III) as one `E8` giant spin and computes longitudinal and transverse relaxation times while its axial zero-field-splitting coefficient varies. The scalar Zeeman parameter is 1.9918 and the applied field is 9.40 T. Redfield relaxation uses `tau_c={1e-15}`, identified in the source as the ligand-cage vibrational correlation time; `rlx_keep='labframe'` and `equilibrium='zero'` are set.

The sweep contains ten linearly spaced `B20` values from 0.1 to 10 cm^-1. At each value, the rank-2 giant-spin coefficient array has `icm2hz(B20)` in its B20 position and zeros in the other supplied positions; the other coefficient array is zero. Both giant-spin Euler-angle arrays are zero. The basis is `sphten-liouv` with `approximation='none'`.

## Calculation and plotted output

For each point, the example builds the spin system, evaluates `R=relaxation(spin_system)`, and obtains the `E8` longitudinal and raising states `Lz` and `Lp`. The relaxation times are the negative reciprocal of the corresponding normalised projections: `T1 = -1 / ((Lz' * R * Lz) / (Lz' * Lz))` and `T2 = -1 / ((Lp' * R * Lp) / (Lp' * Lp))`.

It plots both series against B20 with point markers and a logarithmic y-axis. The axes identify zero-field splitting in cm^-1 and relaxation time in seconds. The source estimates a calculation time of minutes; it supplies no tabulated numerical relaxation-time results.
