# examples/esr_liq_pulsed/pulse_acquire_phenyl.m

## Purpose

A W-band pulse-acquire FFT ESR example for the phenyl radical. It takes the phenyl spin-system data from a vacuum-DFT log and uses simple fixed-linewidth diagonal damping for relaxation.

## Interface and input

- Call from MATLAB with no arguments: `pulse_acquire_phenyl()`. The function declares no return values; it displays a figure rather than saving or returning the FID or spectrum.
- The input path is `../standard_systems/phenyl.log`, parsed by `gparse` and converted by `g2spinach` with electron and proton isotope mappings. `options.no_xyz=1` tells the importer to ignore coordinates because the hyperfine couplings are supplied in the input.

## Spin system and experiment

The source sets `sys.magnet=3.5` and uses the full `sphten-liouv` basis (`approximation='none'`), with the longitudinal `1H` component and projection `+1`. Relaxation is `damp`, retaining diagonal terms, with zero equilibrium and `inter.damp_rate=1e7` (no rate unit is written in the source).

The ESR acquisition uses `liquid(spin_system,@acquire,parameters,'esr')`: the initial state and receiver are both the electron `L+` operator, with no decoupling and zero offset. It sets sweep `2e8`, 512 acquired points, 1024-point zero filling, GHz-labframe axis units, derivative display, and inverted axis. No apodisation is applied; the FID is Fourier transformed and the real spectrum is plotted.

## Dependencies and caveats

Requires Spinach import, system/basis, ESR acquisition, apodisation, and plotting functions (`gparse`, `g2spinach`, `create`, `basis`, `state`, `liquid`, `acquire`, `apodisation`, `kfigure`, `plot_1d`). Run where the relative phenyl-log path resolves. The source supplies no numeric hyperfine table; those couplings come from the log. The magnetic-field value is assigned as 3.5 without an inline unit annotation.

[Source: `pulse_acquire_phenyl.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_liq_pulsed/pulse_acquire_phenyl.m).
