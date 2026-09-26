# examples/dnp_sol/steady_state/tppm_q_con_time_single.m

- Signature: `tppm_q_con_time_single()`

## Purpose

Simulation of TPPM DNP contact time dependence in the steady state. Calculation time: minutes.

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Simulation of TPPM DNP contact time dependence in the
- steady state.
- Calculation time: minutes.
- Q-band magnet
- Electron and proton
- Zeeman interactions (g-tensor for trityl, ppm guess for 1H)
- Spin temperature
- Cartesian coordinates
- Get electron-nuclear distance
- Relaxation rates, distance and ori. dep. R1n
- Basis set
- Propagator accuracy
