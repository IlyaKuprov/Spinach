# examples/giant_spin/quartet_magn.m

- MATLAB implementation: [examples/giant_spin/quartet_magn.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/giant_spin/quartet_magn.m)

- Signature: `quartet_magn()`
- Source: `examples/giant_spin/quartet_magn.m`

## Model and parameters

This example calculates sample magnetisation during a finite-speed field sweep for a spin-3/2 particle (`E4`) with zero-field splitting. It sets `sys.magnet=1.0` T and uses the isotropic Zeeman matrix `diag([2 2 2])`. The splitting parameters are `D=icm2hz(-0.5)` and `E=0.3*D`; `zfs2mat(D,E,0,0,0)` constructs the coupling matrix. The `icm2hz` API takes `-0.5` in inverse centimetres (cm⁻¹) and converts the splitting to Hz. Temperature is set to `1.0`; no temperature unit is stated. The basis is `zeeman-hilb` with `approximation='none'`.

The scan covers fields from 0 to 1 with 1000 points over `sweep_time=1e-9` seconds, at orientation `[0 0 0]`, with four states. The function calls `[fields,z_magn]=fieldscan_magn(spin_system,parameters)` and plots `z_magn` against `fields`. The axes label field in Tesla and sample magnetisation. The source estimates a calculation time of seconds; it gives no fixed numerical trace.
