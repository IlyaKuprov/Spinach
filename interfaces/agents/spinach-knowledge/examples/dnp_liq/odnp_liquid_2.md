# examples/dnp_liq/odnp_liquid_2.m

- Signature: `odnp_liquid_2()`
- Calculation time: seconds

## Purpose

Simulates liquid-state Overhauser DNP after a nominal 180-degree electron inversion pulse, then follows the electron and two proton longitudinal signals as they relax and exchange polarisation. No continuous microwave drive is applied during the evolution.

## Spin system and preparation

The model is the same three-spin geometry as `odnp_liquid_1`: two protons and one electron at 3.4 T, with the nuclei at (0, 0, 0) and (0, 2, 0) Å and the electron at (0, 0, 1.5) Å. It uses explicit Zeeman tensors, a complete sphten-liouv basis, Redfield relaxation, Di Bari equilibrium, secular relaxation retention, 298 K, and a 10 ps correlation time. From the thermal-equilibrium state, the script applies `step` with the electron `Lx` operator and an angle of pi radians, and passes the resulting state as `rho0` to the simulation.

## Evolution and output

The ESR-context `liquid` calculation calls `dnp_time_dep` with zero microwave power and offset, a 1 μs time step, and 1000 steps. The plots show the real electron longitudinal signal and the two proton longitudinal signals over 0–1000 μs.
