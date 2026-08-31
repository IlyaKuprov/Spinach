# examples/dnp_sol/steady_state/xix_q_con_time_ensemble_r_T2e.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/steady_state/xix_q_con_time_ensemble_r_T2e.m`
- Signature: `xix_q_con_time_ensemble_r_T2e()`
- Total lines: 131

## Purpose

Simulation of T2e dependence of XiX DNP contact curves in the steady state with electron-proton distance ensemble. Calculation time: hours

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file also defines local helper function(s): `xix_contact_curve_ensemble_r()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Electron relaxation times, seconds; implemented by `T2e=[50e-6 15e-6 5e-6 1.5e-6 0.5e-6]`.
- Lines 16-17: Get the figure started; implemented by `kfigure(); hold on; kgrid`.
- Lines 22-23: Plot the curves; implemented by `for n=1:numel(T2e)`.
- Lines 27-30: Add the legend and save the plot; implemented by `klegend({'$T_{2e}$ = 50 $\mu$s', '$T_{2e}$ = 15 $\mu$s', '$T_{2e}$ = 5 $\mu$s','$T_{2e}$ = 1.5 $\mu$s', '$T_{2e}$ = 0.5 $\mu$s'},'Location','Best')`.

### Control flow inferred from the code

- Line 23: `for` loop over `n=1:numel(T2e)`.

### Key state/data transformations

- Lines 14: computes `T2e` using `T2e=[50e-6 15e-6 5e-6 1.5e-6 0.5e-6]`.
- Lines 28-30: computes `klegend({'$T_{2e}$` using `klegend({'$T_{2e}$ = 50 $\mu$s', '$T_{2e}$ = 15 $\mu$s', '$T_{2e}$ = 5 $\mu$s','$T_{2e}$ = 1.5 $\mu$s', '$T_{2e}$ = 0.5 $\mu$s'},'Location','Best')`.
- Lines 29-30: computes `'$T_{2e}$` using `'$T_{2e}$ = 5 $\mu$s','$T_{2e}$ = 1.5 $\mu$s', '$T_{2e}$ = 0.5 $\mu$s'},'Location','Best')`.

### Local helper functions

- Line 36: `xix_contact_curve_ensemble_r()` — `function xix_contact_curve_ensemble_r(T2e)`. Q-band magnet
  - Representative operation: `sys.magnet=1.2142`.
  - Representative operation: `sys.isotopes={'E','1H'}`.

## Implementation structure

- Simulation of T2e dependence of XiX DNP contact
- curves in the steady state with electron-proton
- distance ensemble.
- Calculation time: hours
- Electron relaxation times, seconds
- Get the figure started
- Plot the curves
- Add the legend and save the plot
- Simulation for a specific T2e
- Q-band magnet
- Electron and proton
- Zeeman interactions (g-tensor for trityl, ppm guess for 1H)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `kfigure()`, `kxlabel()`, `time()`, `kylabel()`, `ylim()`, `xix_contact_curve_ensemble_r()`, `T2e()`, `klegend()`, `savefig()`, `gaussleg()`, `r1n_dnp()`, `create()`, `basis()`, `state()`, `loop_counts()`, `dnp()`.
