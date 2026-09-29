# examples/dnp_mas/cross_effect_mas_powder.m

- MATLAB implementation: [examples/dnp_mas/cross_effect_mas_powder.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_mas/cross_effect_mas_powder.m)

- Function: `cross_effect_mas_powder()`

## Purpose

Calculates the steady-state cross-effect DNP enhancement for a powder under MAS with Spinach's `masdnp` routine, following [Mentink-Vigier et al.](https://doi.org/10.1016/j.jmr.2015.07.001). The source notes that Spinach rotation conventions differ from the paper and estimates a runtime of minutes.

## Spin system and relaxation

The model uses `{'E','E','1H'}` and `sys.magnet=9.394`. Both electron Zeeman tensors have principal values `[2.0094 2.0060 2.0017]`; the second has Euler angles `pi*[107 108 124]/180` and the first has zero Euler angles. The proton Zeeman values are zero. The electron–electron coupling is `[23.0e6 -11.5e6 -11.5e6]` with Euler angles `pi*[0 135 0]/180`; electron 1–proton coupling is `[1.5e6 -0.75e6 -0.75e6]` with zero Euler angles.

Nottingham relaxation is specified by `nott_r1e=1/0.3e-3`, `nott_r1n=1/4.0`, `nott_r2e=1/1.0e-6`, and `nott_r2n=1/0.2e-3`; the temperature setting is `100`, equilibrium is `'dibari'`, and retained relaxation terms are `'secular'`. The basis is the complete spherical-tensor Liouville basis (`formalism='sphten-liouv'`, `approximation='none'`).

## Powder calculation and output

The call `masdnp(spin_system,parameters)` uses electron spins `{'E'}`, rotor axis `[sqrt(2/3) 0 sqrt(1/3)]`, rate `12.5e3`, and `max_rank=800`. Microwave settings are `mw_pwr=2*pi*0.85e6`, `mw_frq=-263.366e9`, and `mw_time=1.0`; powder orientations use `'rep_2ang_100pts_sph'`. The proton detection coil is set to `state(spin_system,'Lz','1H')`, and verbosity is zero. The returned value is printed as the steady-state DNP enhancement. The example reports a powder enhancement rather than a plotted rotor-period trajectory.
