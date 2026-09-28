# examples/dnp_sol/solid_effect_timedep_1.m

- Signature: `solid_effect_timedep_1()`

## Purpose

Simulates solid-effect DNP dynamics for an electron coupled to a linear chain of three protons at 7, 10, and 14 Å. The chain is tilted by the rotation generated from Euler angles [pi/6, pi/7, pi/8]. The source reports a runtime of seconds.

## Model and method

The calculation uses a 3.4 T field and Weizmann relaxation with secular retention, IME equilibrium, temperature 4.2, and explicitly specified electron/nuclear relaxation rates. Distance-dependent rates are zero except for symmetric 0.1 entries between adjacent proton pairs 1–2 and 2–3. The basis is `sphten-liouv`, with no approximation and projections [-2, -1, 0, 1, 2].

A 250 kHz microwave power and 144.76 MHz nuclear frequency are used. The solid-effect calculation selects `kb_second_order` theory and uses 0.01 s time steps for 1000 steps. It first requests time dependence, then requests a steady state.

## Outputs

The time-domain plots show the real longitudinal expectation value of the electron and of each proton, with logarithmic time axes. The steady-state result is then printed as real `Tr(Sz*rho)` values for all spins.
