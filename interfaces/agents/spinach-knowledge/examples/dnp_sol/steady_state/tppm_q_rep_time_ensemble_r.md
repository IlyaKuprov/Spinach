# examples/dnp_sol/steady_state/tppm_q_rep_time_ensemble_r.m

- Signature: `tppm_q_rep_time_ensemble_r()`

## Purpose

Simulation of TPPM DNP repetition time scan in the steady state with distributions in electron-proton distance. Calculation time: minutes.

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- Simulation of TPPM DNP repetition time scan in the steady
- state with distributions in electron-proton distance.
- Calculation time: minutes.
- Q-band magnet
- Electron and proton
- Zeeman interactions (g-tensor for trityl, ppm guess for 1H)
- Spin temperature
- Basis set
- Propagator accuracy
- Algorithmic options
- Distance ensemble
- Log spacing for rep. time
