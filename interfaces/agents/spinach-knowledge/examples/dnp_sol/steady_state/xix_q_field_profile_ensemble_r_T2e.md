# examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_r_T2e.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_r_T2e.m`
- Signature: `xix_q_field_profile_ensemble_r_T2e()`
- Total lines: 116

## Purpose

Simulation of T2e dependence of XiX DNP field profiles in the steady state with electron-proton distance ensemble. Calculation time: minutes

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file also defines local helper function(s): `xix_field_profile_ensemble_r()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Implementation structure

- Simulation of T2e dependence of XiX DNP field
- profiles in the steady state with electron-proton
- distance ensemble.
- Calculation time: minutes
- Electron relaxation times, seconds
- Get the figure started
- Plot the curves
- Add the legend and save the plot
- Simulation for a specific T2e
- Q-band magnet
- Electron and proton
- Zeeman interactions (g-tensor for trityl, ppm guess for 1H)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `kfigure()`, `kxlabel()`, `kylabel()`, `ylim()`, `xix_field_profile_ensemble_r()`, `T2e()`, `klegend()`, `savefig()`, `gaussleg()`, `r1n_dnp()`, `create()`, `basis()`, `state()`, `dnp()`, `powder()`.
