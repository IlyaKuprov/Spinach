# examples/dnp_mas/solid_effect_mas_dynam.m

- Signature: `solid_effect_mas_dynam()`

## Purpose

Tracks the level populations during one rotor period for a single-crystal solid-effect DNP model. The example follows Fred Mentink-Vigier et al.; Spinach uses different rotation conventions ([paper](https://doi.org/10.1016/j.jmr.2015.07.001)). The source estimates seconds to run.

## Model

The model is an electron–`^1H` pair at 9.403 T. It sets the electron g-tensor principal values to [2.00614, 2.00194, 2.00988] with Euler angles `pi*[253.6 105.1 123.8]/180`, and places the spins 3.00 Å apart on z. Weizmann relaxation parameters are used at 100 K with DiBari equilibrium and secular relaxation retention.

## Calculation

The script builds a full sphten-Liouville basis and an ESR rotor stack about `[sqrt(2/3) 0 sqrt(1/3)]`, with 12.5 kHz MAS and rank limit 3000. It constructs electron `Lx` microwave and `Lz` offset operators, then steps the equilibrium state through the rotor stack with 0.85 MHz microwave power and −400 MHz offset. Finally, `trajan(...,'level_populations')` analyses the resulting one-period trajectory. This is a single-crystal trajectory calculation, not a powder average.
