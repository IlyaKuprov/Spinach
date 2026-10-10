# examples/dnp_mas/solid_effect_mas_dynam.m

- MATLAB implementation: [examples/dnp_mas/solid_effect_mas_dynam.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_mas/solid_effect_mas_dynam.m)

- Function: `solid_effect_mas_dynam()`

## Purpose

Tracks spin-level populations during the first rotor period for a single-crystal solid-effect MAS DNP model. The example follows [Mentink-Vigier et al.](https://doi.org/10.1016/j.jmr.2015.07.001); its source notes that Spinach rotation conventions differ from the paper and estimates seconds of calculation time.

## Spin system and relaxation

The spins are `{'E','1H'}`, with `sys.magnet=9.403`. The electron Zeeman principal values are `[2.00614 2.00194 2.00988]`, with Euler angles `pi*[253.6 105.1 123.8]/180`; the proton Zeeman values are zero. The relative coordinates are `[0 0 0]` and `[0 0 3.00]`.

Weizmann relaxation uses `weiz_r1e=1/0.3e-3`, `weiz_r1n=1/4.0`, `weiz_r2e=1/1.0e-6`, and `weiz_r2n=1/0.2e-3`; both `weiz_r1d` and `weiz_r2d` are zero 2-by-2 arrays. Temperature is `100`, equilibrium is `'dibari'`, and retained relaxation terms are `'secular'`. The complete spherical-tensor Liouville basis is used (`formalism='sphten-liouv'`, `approximation='none'`).

## Rotor-period trajectory

The ESR rotor stack uses axis `[sqrt(2/3) 0 sqrt(1/3)]`, orientation `[0 0 0]`, selected spins `{'E','1H'}`, magnetic MAS frame, empty rotor frames, zero offsets, and `max_rank=3000`; the rate setting is `12.5e3`. Electron `Lx` and `Lz` operators provide the microwave and offset terms, with `mw_pwr=0.85e6` and `mw_off=-400e6`; the relaxation superoperator is included during propagation.

The code forms `rho_eq` from the thermal equilibrium of the left lab-frame Hamiltonian and assigns it as the initial state `rho(:,1)`. The code steps it through the rotor stack with step duration `1/(nsteps*parameters.rate)`, where `nsteps=numel(H)`, then uses `trajan(spin_system,rho,'level_populations')` to analyse the trajectory. This example covers one initial-state rotor period; it does not solve for a periodic steady state or report a DNP enhancement, and it is not a powder average.
