# examples/dnp_mas/solid_effect_mas_powder.m

- Signature: `solid_effect_mas_powder()`

## Purpose

Runs a steady-state MAS DNP powder calculation based on Fred Mentink-Vigier et al. (Spinach rotation conventions differ; [paper](https://doi.org/10.1016/j.jmr.2015.07.001)). The source estimates minutes to run.

## Model and calculation

The model is an electron–`^1H` pair at 9.403 T with the specified anisotropic electron g tensor and a 3.00 Å separation. It uses Weizmann relaxation, DiBari equilibrium, secular relaxation retention, and 100 K temperature. The full sphten-Liouville basis is used. The script calls `masdnp` for the electron-spin powder calculation with 12.5 kHz MAS, rank limit 800, microwave power `2*pi*0.85e6`, frequency `-263.366e9`, and the `rep_2ang_100pts_sph` grid. It reports the returned steady-state DNP enhancement.
