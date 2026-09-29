# examples/giant_spin/quartet_levels.m

- MATLAB implementation: [examples/giant_spin/quartet_levels.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/giant_spin/quartet_levels.m)

- Signature: `quartet_levels()`
- Source: `examples/giant_spin/quartet_levels.m`

## Model and parameters

This field scan models a spin-3/2 particle (`E4`) with a zero-field splitting. The example sets the field to 1.0 T and uses the isotropic Zeeman matrix `diag([2 2 2])`. It defines `D=icm2hz(-0.5)` and `E=0.3*D`, then constructs the zero-field-splitting matrix with `zfs2mat(D,E,0,0,0)`. The `icm2hz` API takes `-0.5` in inverse centimetres (cm⁻¹) and converts the splitting to Hz.

The Hilbert-space basis uses `zeeman-hilb` with `approximation='none'`. The scan parameters are fields from 0 to 1, 100 points, orientation `[0 0 0]`, and four states.

## Calculation

The function creates and bases the Spinach system, then calls `fieldscan_enlev(spin_system,parameters)` for the energy-level field scan. The source estimates a calculation time of seconds; it does not assign the function result or specify a separate plot or fixed numerical output.
