# examples/nmr_solids/mas_powder_nqi_gridfree.m

Source: [examples/nmr_solids/mas_powder_nqi_gridfree.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/mas_powder_nqi_gridfree.m)

## Purpose

Calculates a powder magic-angle-spinning NMR spectrum for one quadrupolar deuterium nucleus. The source describes this as grid-free Fokker-Planck MAS and notes that second-order corrections to the rotating-frame transformation are not applied. Its “minutes” runtime is a source comment estimate, not a recorded run time.

## Model and acquisition

The model contains `2H` at `sys.magnet=9.4` T, with a quadrupolar tensor specified by principal values `-1000 -2000 3000` rad/s and Euler angles `0 0 0` (the [Spinach NQI convention](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/agents/spinach-knowledge/kernel/conventions/transforms/ham2nqi.md) gives the quadrupolar tensor in rad/s). The interaction values are model inputs, not measured parameters in this example. The basis is `sphten-liouv`, approximation `none`, with projection `+1`.

The rotor-axis vector is `1 1 1` and the MAS rate is `1000` Hz. The acquisition is set to a `20000` Hz sweep, 512 points, zero filling to 4096, zero offset, and ppm axis labelling. No RF field or pulse sequence is specified. The initial state and receiver are both `L+` on `2H`; there is no separately assigned dipolar interaction in this one-spin model.

## Calculation and display

The code calls `gridfree(spin_system,@acquire,parameters,'nmr')`. It exponentially apodises the calculated FID with parameter 6, Fourier transforms it using the 4096-point zero-fill, then plots the real spectrum with `plot_1d`. This is a simulated model spectrum; the source does not provide an experimental trace or a measured comparison.
