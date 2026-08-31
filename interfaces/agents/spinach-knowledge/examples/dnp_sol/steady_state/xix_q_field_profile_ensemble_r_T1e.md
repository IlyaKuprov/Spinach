# examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_r_T1e.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_r_T1e.m`
- Signature: `xix_q_field_profile_ensemble_r_T1e()`
- Total lines: 116

## Purpose

Simulation of T1e dependence of XiX DNP field profiles in the steady state with electron-proton distance ensemble. Calculation time: minutes

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file also defines local helper function(s): `xix_field_profile_ensemble_r()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Electron relaxation times to use, seconds; implemented by `T1e=[10e-3, 3.0e-3, 1.0e-3, 0.3e-3, 0.1e-3]`.
- Lines 15-16: Get the figure started; implemented by `kfigure(); hold on; kgrid`.
- Lines 21-22: Plot the curves; implemented by `for n=1:numel(T1e)`.
- Lines 26-29: Add the legend and save the plot; implemented by `klegend({'$T_{1e}$ = 10 ms', '$T_{1e}$ = 3.0 ms', '$T_{1e}$ = 1.0 ms','$T_{1e}$ = 0.3 ms', '$T_{1e}$ = 0.1 ms'},'Location','Best')`.

### Control flow inferred from the code

- Line 22: `for` loop over `n=1:numel(T1e)`.

### Key state/data transformations

- Lines 13: computes `T1e` using `T1e=[10e-3, 3.0e-3, 1.0e-3, 0.3e-3, 0.1e-3]`.
- Lines 27-29: computes `klegend({'$T_{1e}$` using `klegend({'$T_{1e}$ = 10 ms', '$T_{1e}$ = 3.0 ms', '$T_{1e}$ = 1.0 ms','$T_{1e}$ = 0.3 ms', '$T_{1e}$ = 0.1 ms'},'Location','Best')`.
- Lines 28-29: computes `'$T_{1e}$` using `'$T_{1e}$ = 1.0 ms','$T_{1e}$ = 0.3 ms', '$T_{1e}$ = 0.1 ms'},'Location','Best')`.

### Local helper functions

- Line 35: `xix_field_profile_ensemble_r()` — `function xix_field_profile_ensemble_r(T1e)`. Q-band magnet
  - Representative operation: `sys.magnet=1.2142`.
  - Representative operation: `sys.isotopes={'E','1H'}`.

## Implementation structure

- Simulation of T1e dependence of XiX DNP field profiles
- in the steady state with electron-proton distance ensemble.
- Calculation time: minutes
- Electron relaxation times to use, seconds
- Get the figure started
- Plot the curves
- Add the legend and save the plot
- Simulation for a specific T1e
- Q-band magnet
- Electron and proton
- Zeeman interactions (g-tensor for trityl, ppm guess for 1H)
- Spin temperature

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `kfigure()`, `kxlabel()`, `kylabel()`, `ylim()`, `xix_field_profile_ensemble_r()`, `T1e()`, `klegend()`, `savefig()`, `gaussleg()`, `r1n_dnp()`, `create()`, `basis()`, `state()`, `dnp()`, `powder()`.
