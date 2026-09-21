# examples/dnp_sol/steady_state/xix_q_nutation_ensemble_b1_r.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/steady_state/xix_q_nutation_ensemble_b1_r.m`
- Signature: `xix_q_nutation_ensemble_b1_r()`
- Total lines: 134

## Purpose

Simulation of nutation frequency dependence of XiX DNP field profiles in the steady state with electron-proton distance and electron Rabi frequency ensembles. Calculation time: minutes

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file also defines local helper function(s): `xix_field_profile_b1_r()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Implementation structure

- Simulation of nutation frequency dependence of XiX DNP
- field profiles in the steady state with electron-proton
- distance and electron Rabi frequency ensembles.
- Calculation time: minutes
- Nutation frequencies, Hz
- Shot repetition times, seconds
- Get the figure started
- Plot the curves
- Save results
- Q-band magnet
- Electron and proton
- Zeeman interactions (g-tensor for trityl, ppm guess for 1H)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `kfigure()`, `kxlabel()`, `kylabel()`, `kzlabel()`, `view()`, `xlim()`, `ylim()`, `zlim()`, `set()`, `xix_field_profile_b1_r()`, `srt()`, `savefig()`, `gaussleg()`, `r1n_dnp()`, `create()`, `basis()`.
