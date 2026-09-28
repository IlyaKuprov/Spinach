# examples/optimal_control/steady_orbit/solid_effect_int.m

- Signature: `solid_effect_int()`

## Purpose

Panoramic optimisation for stroboscopic steady-state DNP, using timing and power settings matching the XiX experiment while allowing the phase to vary freely.

## Model and numerical setup

- An electron (`E`) and proton (`1H`) at a 3.35316 T W-band field (HIPER at St Andrews), separated by 3.500 Å, with a spin temperature of 80 K.
- Zeeman principal values are `[2.00319 2.00319 2.00258]` for the trityl electron and `[0 0 5]` ppm for the proton; Euler angles are `[0 10 0]` and `[0 0 10]` degrees.
- `t1_t2` relaxation uses electron R1 = `1e3`, R2 = `200e3`, proton R2 = `50e3`, and an orientation-dependent proton R1 from `r1n_dnp` using `2.00230`, `1.0e-3`, `52.0`, the electron–proton distance, and `bet`. Relaxation is diagonal and equilibrium is `dibari`.
- The calculation uses an unrestricted `sphten-liouv` basis, 240 processes, propagation chopping tolerance `1e-14`, and steady-state tolerance `1e-10`. Calculation time is days on a large parallel cluster.

## Optimisation

- Powder averaging uses `rep_2ang_800pts_sph`; the transmitter is set to precisely 94.0 GHz. The target is proton `Lz` magnetisation normalised to its thermal-equilibrium value.
- Electron `Lx` and `Ly` are phase controls, with `Lz` as the offset operator. Microwave power levels are `2*pi*linspace(5,25,20)*1e6` rad/s, and tested offsets are `[-2 -1 0 +1 +2]` MHz.
- The sequence comprises 720 adjustable 0.5 ns pulse samples, 20 frozen 0.5 ns ringdown samples, and a frozen 167 µs delay. Amplitude is one during the pulse and zero thereafter.
- The HiPER filter is loaded from `hiper_kernel_trans.mat`; its first 16 coefficients are normalised to unit absolute DC gain and applied through `firf` in optimisation and plotting.
- A sinusoidal phase chirp initialises `fmaxnewton` with `grape_phase`. Optimisation uses `rbfgs`, up to 10,000 iterations, `steady=true`, and a budget of 500; robustness and spectrogram plots are requested.
