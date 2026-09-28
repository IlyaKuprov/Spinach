# examples/dnp_sol/steady_state/xix_q_con_time_single.m

- Signature: `xix_q_con_time_single()`

## Purpose

Simulate the steady-state XiX DNP proton signal as a function of total contact time at a Q-band magnetic field. The source estimates a calculation time of minutes.

## Physical / mathematical content

- Models one electron and one proton (`'E'`, `'1H'`) at `sys.magnet=1.2142` and `inter.temperature=80`. The electron has a trityl g-tensor; the proton has an illustrative Zeeman shift. Their coordinates place them 3.500 units apart along z.
- Uses `t1_t2` relaxation, including a proton longitudinal rate supplied by `r1n_dnp` as a function of orientation angle `bet` and electron–nuclear distance. Equilibrium is set to `'dibari'`, and the detected observable is proton `Lz`.
- Averages the steady-state XiX response over the specified orientation grid using `powder(spin_system,@xixdnp_steady,parameters,'esr')`.

## Numerical / algorithmic content

- Creates the Spinach system with `create` and `basis`, using an unrestricted `sphten-liouv` basis and propagator chopping tolerance `1e-12`.
- Sweeps `parameters.nloops` from 1 to 64. Each loop contributes two 48 ns pulses, so total contact time is `2*pulse_dur*nloops`; shot spacing is updated to `153e-6` minus that pulse-train duration.
- Sets electron irradiation to 18 MHz, second-pulse phase to `pi`, and additional shift and electron offset to −13 MHz and 61 MHz. Plots the real part of the calculated proton `Lz` signal against contact time in microseconds and saves `xix_q_con_time_single.fig`.
