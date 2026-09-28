# examples/giant_spin/lanthanide_redfield.m

- Signature: `lanthanide_redfield()`

## Purpose

Calculate longitudinal and transverse relaxation times for an `E8` giant spin as the zero-field-splitting coefficient $B_2^0$ varies. The calculation uses Redfield relaxation with a 1 fs correlation time, identified in the source as the timescale of ligand-cage vibrations. Estimated calculation time: minutes.

## Physical / mathematical content

- The spin system contains one `E8` spin with scalar Zeeman parameter 1.9918 in a 9.40 T magnetic field.
- For each $B_2^0$, the giant-spin coefficient arrays contain zero entries except for the $B_2^0$ entry, which is converted from cm⁻¹ to Hz with `icm2hz`. Both giant-spin Euler-angle arrays are zero.
- Redfield relaxation is configured with `rlx_keep='labframe'`, `equilibrium='zero'`, and `tau_c=1e-15` s.
- The relaxation times are obtained from the relaxation superoperator `R` using normalized operator projections: $T_1=-1/[(L_z^\dagger R L_z)/(L_z^\dagger L_z)]$ and $T_2=-1/[(L_+^\dagger R L_+)/(L_+^\dagger L_+)]$, where `Lz` and `L+` are states of `E8`.

## Numerical / algorithmic content

- Sweep ten linearly spaced $B_2^0$ values from 0.1 to 10 cm⁻¹.
- At each value, create the spin system, construct a `sphten-liouv` basis with `approximation='none'`, calculate `R=relaxation(spin_system)`, and evaluate $T_1$ and $T_2$.
- Plot both relaxation-time series against $B_2^0$ using point markers and a logarithmic vertical axis. The horizontal axis is labeled in cm⁻¹ and the vertical axis in seconds.

## Implementation structure

`lanthanide_redfield()` sets the spin, field, relaxation, and basis parameters; loops over the $B_2^0$ sweep; computes the two relaxation times; and plots the results. It takes no arguments and returns no values.
