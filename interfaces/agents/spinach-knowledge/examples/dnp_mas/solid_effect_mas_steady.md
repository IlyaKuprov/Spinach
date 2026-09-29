# examples/dnp_mas/solid_effect_mas_steady.m

- MATLAB implementation: [examples/dnp_mas/solid_effect_mas_steady.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_mas/solid_effect_mas_steady.m)

## Purpose

Call `solid_effect_mas_steady()` to find the steady state over one MAS rotor period for the single-crystal electron–`^1H` solid-effect DNP model, then inspect the rotor-period population trajectory. The example follows Fred Mentink-Vigier's paper (Spinach's rotation conventions differ; [paper](https://doi.org/10.1016/j.jmr.2015.07.001)); its header estimates seconds.

## Model and calculation

The system has field `9.403`, electron g eigenvalues `[2.00614 2.00194 2.00988]` and Euler angles `pi*[253.6 105.1 123.8]/180`, zero proton Zeeman tensor, and coordinates `[0 0 0]` and `[0 0 3.00]`. Weizmann relaxation is set to `weiz_r1e=1/0.3e-3`, `weiz_r1n=1/4.0`, `weiz_r2e=1/1.0e-6`, and `weiz_r2n=1/0.2e-3`, with zero dipolar relaxation matrices. Temperature is `100`, equilibrium is `'dibari'`, and relaxation retention is `'secular'`. The basis is the full sphten-Liouville basis.

The ESR rotor stack uses axis `[sqrt(2/3) 0 sqrt(1/3)]`, empty rotor frames, orientation `[0 0 0]`, spins `{'E','1H'}`, magnet MAS frame, zero offsets, and `max_rank=3000`. The microwave terms use `parameters.mw_pwr=0.85e6` and `parameters.mw_off=-400e6`; the source combines them with `2*pi` in the Hamiltonian. It creates `Hmw=operator(spin_system,'Lx','E')` and `HzE=operator(spin_system,'Lz','E')` and obtains the relaxation superoperator `R=relaxation(spin_system)`. For each stack element, it composes the rotor-period propagator using timestep `1/(nsteps*parameters.rate)`, then calls `steady(spin_system,P,[],'newton')`.

## Output and scope

Starting from that steady state, the code steps through the rotor-stack elements and calls `trajan(spin_system,rho,'level_populations')`. This yields the level-population trajectory for analysis. The script does not call powder averaging or print a separate scalar proton enhancement.
