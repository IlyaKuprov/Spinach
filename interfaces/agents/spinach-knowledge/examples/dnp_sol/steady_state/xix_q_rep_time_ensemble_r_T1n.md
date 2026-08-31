# examples/dnp_sol/steady_state/xix_q_rep_time_ensemble_r_T1n.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/steady_state/xix_q_rep_time_ensemble_r_T1n.m`
- Signature: `xix_q_rep_time_ensemble_r_T1n()`
- Total lines: 126

## Purpose

Simulation of T1n dependent XiX DNP optimisation of repetition time in the steady state with electron- proton distance ensemble. Calculation time: hours

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file also defines local helper function(s): `xix_rep_time_ensemble_r()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Nuclear relaxation times, seconds; implemented by `T1n=[50.0 5.0 0.500 0.050 0.005]`.
- Lines 16-17: Get the figure started; implemented by `kfigure(); hold on; kgrid`.
- Lines 22-23: Plot the curves; implemented by `for n=1:numel(T1n)`.
- Lines 27-30: Add the legend and save the plot; implemented by `klegend({'$T_{1n}$ = 50 s', '$T_{1n}$ = 5 s', '$T_{1n}$ = 0.5 s','$T_{1n}$ = 0.05 s', '$T_{1n}$ = 0.005 s'},'Location','Best')`.

### Control flow inferred from the code

- Line 23: `for` loop over `n=1:numel(T1n)`.

### Key state/data transformations

- Lines 14: computes `T1n` using `T1n=[50.0 5.0 0.500 0.050 0.005]`.
- Lines 28-30: computes `klegend({'$T_{1n}$` using `klegend({'$T_{1n}$ = 50 s', '$T_{1n}$ = 5 s', '$T_{1n}$ = 0.5 s','$T_{1n}$ = 0.05 s', '$T_{1n}$ = 0.005 s'},'Location','Best')`.
- Lines 29-30: computes `'$T_{1n}$` using `'$T_{1n}$ = 0.5 s','$T_{1n}$ = 0.05 s', '$T_{1n}$ = 0.005 s'},'Location','Best')`.

### Local helper functions

- Line 36: `xix_rep_time_ensemble_r()` — `function xix_rep_time_ensemble_r(T1n)`. Q-band magnet
  - Representative operation: `sys.magnet=1.2142`.
  - Representative operation: `sys.isotopes={'E','1H'}`.

## Implementation structure

- Simulation of T1n dependent XiX DNP optimisation of
- repetition time in the steady state with electron-
- proton distance ensemble.
- Calculation time: hours
- Nuclear relaxation times, seconds
- Get the figure started
- Plot the curves
- Add the legend and save the plot
- Simulation for a specific T1n
- Q-band magnet
- Electron and proton
- Zeeman interactions (g-tensor for trityl, ppm guess for 1H)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `kfigure()`, `kxlabel()`, `time()`, `kylabel()`, `xlim()`, `ylim()`, `xix_rep_time_ensemble_r()`, `T1n()`, `klegend()`, `savefig()`, `gaussleg()`, `create()`, `basis()`, `state()`, `rep_time()`, `dnp()`.
