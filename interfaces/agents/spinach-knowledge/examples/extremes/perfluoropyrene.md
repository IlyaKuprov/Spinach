# examples/extremes/perfluoropyrene.m

- MATLAB implementation: [examples/extremes/perfluoropyrene.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/extremes/perfluoropyrene.m)

- Entry point: `perfluoropyrene()`.
- Source: [`examples/extremes/perfluoropyrene.m`](../../../../../examples/extremes/perfluoropyrene.m); input log: [`examples/standard_systems/perfluoropyrene_cation.log`](../../../../../examples/standard_systems/perfluoropyrene_cation.log).

## Intent and spin-system construction

The source describes this as an X-band pulsed ESR spectrum of the perfluoropyrene cation radical. It reads the molecular data from the linked log using `gparse` and `g2spinach`, with the supplied electron/`19F` isotope mapping; `options.no_xyz=1` tells the conversion to ignore coordinate information. The explicit input and conversion settings, rather than a hand-written isotope list, define the spin system. The magnetic-induction parameter is `0.33` (unit not written in the script).

Relaxation is configured as `{'damp'}`, with diagonal relaxation operators, zero equilibrium state and `damp_rate=2e6`. The basis is `sphten-liouv` with `approximation='none'`. The source comment calls this brute-force operator algebra in the full 4,194,304-dimensional Liouville space and notes that a restricted basis would be faster; it also identifies trajectory-level state-space restriction as the performance technique of interest.

## ESR acquisition and output

The observed spin is `E` (electron), with `L+` initial state and receiver, no decoupled spins, and offset `0`. The sweep is `3e8`, with `2048` points and zero-fill to `4096`; the axis is `GHz-labframe`, the derivative parameter is `1`, and the axis is inverted. The script does not give units for offset, sweep or damping rate beyond the axis label.

It runs `liquid(spin_system,@acquire,parameters,'esr')`, applies no apodisation (`{'none'}`), Fourier-transforms the FID and plots the real spectrum. It defines no explicit pulse-duration or pulse-amplitude table; the sequence is the acquisition callback invoked by the ESR-mode engine.

## Practical limit

The source comments state a minimum of 64 GB RAM and a calculation time of minutes.
