# examples/dnp_mas/solid_effect_mas_enlev.m

- Signature: `solid_effect_mas_enlev()`

## Purpose

Plots the energy levels of a single-crystal electron–`^1H` model against MAS rotor phase, following Fred Mentink-Vigier et al. (Spinach rotation conventions differ; [paper](https://doi.org/10.1016/j.jmr.2015.07.001)). The source estimates milliseconds to run.

## Model and calculation

At 9.403 T, the example defines the electron g tensor and the 3.00 Å spin separation, then uses the full Zeeman-Hilbert basis. It generates an ESR rotor stack about `[sqrt(2/3) 0 sqrt(1/3)]` with the magnet as MAS frame and rank limit 200. For each stack Hamiltonian, a `parfor` loop sorts the real parts of its eigenvalues. The results are plotted against rotor phase in radians, with energy levels in rad/s. This is an energy-level plot, not a powder-averaged DNP calculation.
