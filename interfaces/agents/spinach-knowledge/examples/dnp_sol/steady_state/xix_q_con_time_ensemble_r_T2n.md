# examples/dnp_sol/steady_state/xix_q_con_time_ensemble_r_T2n.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/steady_state/xix_q_con_time_ensemble_r_T2n.m`
- Signature: `xix_q_con_time_ensemble_r_T2n()`
- Total lines: 131

## Purpose

Simulation of T2n dependence of XiX DNP contact curves in the steady state with electron-proton distance ensemble. Calculation time: hours

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file also defines local helper function(s): `xix_contact_curve_ensemble_r()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Implementation structure

- Simulation of T2n dependence of XiX DNP contact
- curves in the steady state with electron-proton
- distance ensemble.
- Calculation time: hours
- Nuclear T2 times, seconds
- Get the figure started
- Plot the curves
- Add the legend and save the plot
- Simulation for a specific T2n
- Q-band magnet
- Electron and proton
- Zeeman interactions (g-tensor for trityl, ppm guess for 1H)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `kfigure()`, `kxlabel()`, `time()`, `kylabel()`, `ylim()`, `xix_contact_curve_ensemble_r()`, `T2n()`, `klegend()`, `savefig()`, `gaussleg()`, `r1n_dnp()`, `create()`, `basis()`, `state()`, `loop_counts()`, `dnp()`.
