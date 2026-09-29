# examples/dnp_mas/cross_effect_mas_steady.m

- MATLAB implementation: [examples/dnp_mas/cross_effect_mas_steady.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_mas/cross_effect_mas_steady.m)

- Function: `cross_effect_mas_steady()`

## Purpose

Finds the periodic steady state of a single-crystal cross-effect MAS DNP model, analyses its trajectory over a rotor period, and reports the proton enhancement. It follows [Mentink-Vigier et al.](https://doi.org/10.1016/j.jmr.2015.07.001); the source notes that Spinach rotation conventions differ from the paper and estimates minutes of calculation time.

## Spin system and rotor stack

The spins are `{'E','E','1H'}`, with `sys.magnet=9.394`. Both electron Zeeman tensors have principal values `[2.0094 2.0060 2.0017]`; the second uses Euler angles `pi*[107 108 124]/180`, while the first uses zero angles. The proton Zeeman values are zero. Electron–electron coupling is `[23.0e6 -11.5e6 -11.5e6]` with Euler angles `pi*[0 135 0]/180`; electron 1–proton coupling is `[1.5e6 -0.75e6 -0.75e6]` with zero angles.

The complete spherical-tensor Liouville basis is used (`formalism='sphten-liouv'`, `approximation='none'`). Nottingham relaxation has settings `nott_r1e=1/0.3e-3`, `nott_r1n=1/4.0`, `nott_r2e=1/1.0e-6`, and `nott_r2n=1/0.2e-3`; temperature is `100`, equilibrium is `'dibari'`, and relaxation retention is `'secular'`.

The ESR rotor stack uses axis `[sqrt(2/3) 0 sqrt(1/3)]`, orientation `pi*[320 141 80]/180`, selected spins `{'E','1H'}`, magnetic MAS frame, empty rotor frames, zero offsets, and `max_rank=3000`. The rotor rate setting is `12.5e3`. Electron microwave drive and offset settings are `mw_pwr=0.85e6` and `mw_off=-400e6`, applied through electron `Lx` and `Lz` operators; the relaxation superoperator is included in each interval generator.

## Periodic state, trajectory, and observable

The code composes one-interval propagators around the rotor stack and calls `steady(spin_system,P,[],'newton')` for the periodic state. It then steps that state through the rotor intervals and plots `trajan(spin_system,rho,'level_populations')`. The reported enhancement is the real ratio of the period-mean proton `Lz` signal to the proton `Lz` signal in a thermal-equilibrium state formed from the lab-frame Hamiltonian. This is a fixed-orientation single-crystal calculation, not a powder average.
