# examples/esr_liq_pulsed/relaxation_bisnitroxide.m

## Purpose and interface

An X-band pulse-acquire FFT ESR example for a bisnitroxide, using explicit time-domain simulation with a Redfield relaxation superoperator. Call `relaxation_bisnitroxide()` with no arguments. It declares no return values and plots the real spectrum; it does not save the signal or spectrum.

## Spin system and model

The spins are two electrons and two `14N` nuclei (`{'E','E','14N','14N'}`). Both electron principal g values are `[2.00925 2.00605 2.00205]`; the first orientation is `[0 0 0]` and the second is `[123.1 129.8 -46.6]*pi/180`. The nitrogen Zeeman tensors are set to zero. Coupling eigenvalues are encoded as electron-electron `[17.5 17.5 -35]*1e6` and electron-nitrogen `[18 17 103]*1e6` (for pairs 1–2 and 1–3 / 2–4 respectively); the source does not state their units. The 1–2 tensor Euler angles are `[-174 74 0]*pi/180`, the 1–3 angles are zero, and the 2–4 angles copy electron 2's Zeeman orientation. The scalar 1–2 coupling is assigned `-2*16e6`, also without a source unit label.

The source sets `sys.magnet=0.35` and comments this value as Tesla. It uses the full `sphten-liouv` basis and Redfield relaxation with zero equilibrium, lab-frame terms retained, and `inter.tau_c={4e-10}` (no unit is attached to this value in the source).

## Pulse-acquire and processing

Initial density and receiver are both electron `L+`; decoupling is empty and offset is `0e8`. Acquisition uses sweep `5e8`, 512 points, 1024-point zero filling, GHz-labframe axis units, first-derivative display, and inverted axis. The `liquid(...,@acquire,...,'esr')` signal is not apodised, then Fourier transformed and plotted. Sweep/offset and coupling values are left in their source-coded form because no units are given for them.

Requires Spinach `create`, `basis`, `state`, `liquid`, `acquire`, `apodisation`, `kfigure`, and `plot_1d`. The parameter source cited by the example is [DOI: 10.1039/C8CP06819D](https://doi.org/10.1039/C8CP06819D).

[Source: `relaxation_bisnitroxide.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_liq_pulsed/relaxation_bisnitroxide.m).
