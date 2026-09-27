# examples/dnp_mas/solid_effect_mas_steady.m

- Signature: `solid_effect_mas_steady()`

## Purpose

Finds the steady state of a single-crystal solid-effect DNP system over one MAS rotor period, then analyses its level-population trajectory and proton enhancement. The example follows Fred Mentink-Vigier et al. (Spinach rotation conventions differ; [paper](https://doi.org/10.1016/j.jmr.2015.07.001)); the source estimates seconds to run.

## Model and calculation

The model is an electron–`^1H` pair at 9.403 T with the specified anisotropic electron g tensor and 3.00 Å separation. It uses Weizmann relaxation at 100 K, DiBari equilibrium, secular relaxation retention, and a full sphten-Liouville basis. An ESR rotor stack (12.5 kHz MAS, rank limit 3000) is formed about `[sqrt(2/3) 0 sqrt(1/3)]`. With 0.85 MHz microwave power and −400 MHz offset, the code composes one-period propagator `P`, finds `rho_st=steady(...,'newton')`, steps that state through the rotor period, and calls `trajan` for level populations. It reports the proton `Lz` expectation relative to thermal equilibrium as the enhancement factor. This is a single-crystal calculation, not a powder average.
