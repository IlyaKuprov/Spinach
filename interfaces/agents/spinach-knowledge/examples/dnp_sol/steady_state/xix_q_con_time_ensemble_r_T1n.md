# examples/dnp_sol/steady_state/xix_q_con_time_ensemble_r_T1n.m

- Signature: `xix_q_con_time_ensemble_r_T1n()`

## Purpose

Simulation of T1n dependence of XiX DNP contact curves in the steady state with electron-proton distance ensemble. Calculation time: hours

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Simulation of T1n dependence of XiX DNP contact
- curves in the steady state with electron-proton
- distance ensemble.
- Calculation time: hours
- Nuclear T1 times, seconds
- Get the figure started
- Plot the curves
- Add the legend and save the plot
- Simulation for a specific T1n
- Q-band magnet
- Electron and proton
- Zeeman interactions (g-tensor for trityl, ppm guess for 1H)
