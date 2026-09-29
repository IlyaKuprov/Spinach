# examples/dnp_mas/solid_effect_mas_powder.m

- MATLAB implementation: [examples/dnp_mas/solid_effect_mas_powder.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_mas/solid_effect_mas_powder.m)

## Purpose

Call the no-argument function `solid_effect_mas_powder()` to calculate and print a steady-state MAS DNP enhancement for a powder. The example follows Fred Mentink-Vigier's treatment; the source notes that Spinach's rotation conventions differ ([paper](https://doi.org/10.1016/j.jmr.2015.07.001)). Its header estimates minutes.

## Model and spin-dynamics settings

The two-spin system is `{'E','1H'}` with field `9.403`, electron g eigenvalues `[2.00614 2.00194 2.00988]` and Euler angles `pi*[253.6 105.1 123.8]/180`; the proton Zeeman tensor is zero. Coordinates are `[0 0 0]` and `[0 0 3.00]`. Relaxation uses the Weizmann model with `weiz_r1e=1/0.3e-3`, `weiz_r1n=1/4.0`, `weiz_r2e=1/1.0e-6`, and `weiz_r2n=1/0.2e-3`; both dipolar relaxation matrices are zero. Temperature is assigned `100`, equilibrium is `'dibari'`, and relaxation retention is `'secular'`. The basis is the full sphten-Liouville basis (`formalism='sphten-liouv'`, `approximation='none'`).

The calculation sets the MAS rate to `12.5e3`, rotor axis to `[sqrt(2/3) 0 sqrt(1/3)]`, electron spins for the stack (`parameters.spins={'E'}`), and rank limit `800`. Microwave parameters are `mw_pwr=2*pi*0.85e6`, `mw_frq=-263.366e9`, and `mw_time=1.0`. Powder averaging uses `'rep_2ang_100pts_sph'`; detection is the proton state `state(spin_system,'Lz','1H')`. The function passes these settings to `masdnp(spin_system,parameters)`.

## Output and scope

The returned scalar is displayed with the label “Steady state DNP enhancement”. This is the powder calculation; unlike the companion energy-level example, it returns an enhancement rather than plotting rotor-phase eigenvalues.
