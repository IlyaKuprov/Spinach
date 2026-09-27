# examples/dnp_liq/odnp_liquid_1.m

- Signature: `odnp_liquid_1()`
- Calculation time: seconds

## Purpose

Demonstrates liquid-state Overhauser DNP for two protons coupled dipolarly to an electron, with continuous-wave electron irradiation on resonance. The time-dependent simulation plots the electron and proton longitudinal signals during irradiation.

## Spin system and experiment

The three spins are two 1H nuclei and one electron at 3.4 T. The protons are at (0, 0, 0) and (0, 2, 0) Å; the electron is at (0, 0, 1.5) Å. Their Zeeman tensors are set explicitly. The model uses a complete sphten-liouv basis, Redfield relaxation, Di Bari equilibrium, secular relaxation retention, temperature 298 K, and a 10 ps correlation time.

The `liquid` simulation requests the equilibrium state, detects all three longitudinal signals, and drives the electron with its Lx operator at zero offset and `parameters.mw_pwr=2*pi*1e6`. With a 1 μs time step and 1000 steps, it calls `dnp_time_dep` in the ESR context and plots the electron signal and the two proton signals over 0–1000 μs.
