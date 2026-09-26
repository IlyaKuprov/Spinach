# examples/dnp_sol/steady_state/top_q_nutation_ensemble_b1_r.m

- Signature: `top_q_nutation_ensemble_b1_r()`

## Purpose

Simulation of nutation frequency dependence of TOP DNP field profiles in the steady state with electron-proton distance and electron Rabi frequency ensembles. Calculation time: minutes

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Simulation of nutation frequency dependence of TOP DNP
- field profiles in the steady state with electron-proton
- distance and electron Rabi frequency ensembles.
- Calculation time: minutes
- Nutation frequencies, Hz
- Get the figure started
- Plot the curves
- Save results
- Q-band magnet
- Electron and proton
- Zeeman interactions (g-tensor for trityl, ppm guess for 1H)
- Spin temperature
