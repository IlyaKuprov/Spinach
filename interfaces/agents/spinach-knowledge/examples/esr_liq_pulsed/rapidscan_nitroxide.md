# examples/esr_liq_pulsed/rapidscan_nitroxide.m

## Purpose and interface

A rapid-scan ESR calculation for a nitroxide radical. Call `rapidscan_nitroxide()` with no arguments in MATLAB. The function has no declared outputs; it plots the calculated real spectrum rather than saving or returning the axis or data.

## Spin system and relaxation

The two-spin system is `{'14N','E'}`, with centre-field setting `sys.magnet=3.5`. Nitrogen's Zeeman matrix is zero; the electron g matrix is

```
[2.0104 0      0.0001
 0      2.0064 0
 0.0001 0      2.0021].
```

The symmetric nitrogen-electron coupling matrix is `[0.6178 0 0.3161; 0 0.5633 0; 0.3161 0 4.1115]*1e7`. The source does not label the units of these tensor entries. The basis is full `sphten-liouv`; relaxation is Redfield with secular terms, `equilibrium='dibari'`, `inter.temperature=100`, and `inter.tau_c={2e-11}` (the source does not attach units to the last two values).

## Rapid-scan settings and output

Unlike the other pages in this group, this one calls `rapidscan(spin_system,parameters)`, not the pulse-acquire `liquid/@acquire` pathway. It sets `mw_pwr=2*pi*1e3`, sweep endpoints `[-0.011 -0.003]`, 500 steps, and timestep `1e-8`. The output is plotted against the returned magnetic-induction axis (labelled T), with signal intensity labelled in arbitrary units; the source does not save the arrays. The power, sweep, and timestep assignments have no unit comments, so they are given here as coded.

Requires Spinach `create`, `basis`, `rapidscan`, and plotting support (`kfigure`, `plot`). No external spin-system file is read.

[Source: `rapidscan_nitroxide.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_liq_pulsed/rapidscan_nitroxide.m).
