# examples/dnp_sol/steady_state/xix_q_rep_time_ensemble_r_T2n.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/steady_state/xix_q_rep_time_ensemble_r_T2n.m`
- Signature: `xix_q_rep_time_ensemble_r_T2n()`
- Total lines: 128

## Purpose

Simulation of T2n dependent XiX DNP optimisation of repetition time in the steady state with electron- proton distance ensemble. Calculation time: hours

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file also defines local helper function(s): `xix_rep_time_ensemble_r()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Nuclear relaxation times, seconds; implemented by `T2n=[2e-3 200e-6 20e-6 2e-6 0.2e-6]`.
- Lines 16-17: Get the figure started; implemented by `kfigure(); hold on; kgrid`.
- Lines 22-23: Plot the curves; implemented by `for n=1:numel(T2n)`.
- Lines 27-30: Add the legend and save the plot; implemented by `klegend({'$T_{2n}$ = 2000 $\mu$s', '$T_{2n}$ = 200 $\mu$s', '$T_{2n}$ = 20 $\mu$s','$T_{2n}$ = 2 $\mu$s', '$T_{2n}$ = 0.2 $\mu$s'},'Location','Best')`.

### Control flow inferred from the code

- Line 23: `for` loop over `n=1:numel(T2n)`.

### Key state/data transformations

- Lines 14: computes `T2n` using `T2n=[2e-3 200e-6 20e-6 2e-6 0.2e-6]`.
- Lines 28-30: computes `klegend({'$T_{2n}$` using `klegend({'$T_{2n}$ = 2000 $\mu$s', '$T_{2n}$ = 200 $\mu$s', '$T_{2n}$ = 20 $\mu$s','$T_{2n}$ = 2 $\mu$s', '$T_{2n}$ = 0.2 $\mu$s'},'Location','Best')`.
- Lines 29-30: computes `'$T_{2n}$` using `'$T_{2n}$ = 20 $\mu$s','$T_{2n}$ = 2 $\mu$s', '$T_{2n}$ = 0.2 $\mu$s'},'Location','Best')`.

### Local helper functions

- Line 36: `xix_rep_time_ensemble_r()` — `function xix_rep_time_ensemble_r(T2n)`. Q-band magnet
  - Representative operation: `sys.magnet=1.2142`.
  - Representative operation: `sys.isotopes={'E','1H'}`.

## Implementation structure

- Simulation of T2n dependent XiX DNP optimisation of
- repetition time in the steady state with electron-
- proton distance ensemble.
- Calculation time: hours
- Nuclear relaxation times, seconds
- Get the figure started
- Plot the curves
- Add the legend and save the plot
- Simulation for a specific T2n
- Q-band magnet
- Electron and proton
- Zeeman interactions (g-tensor for trityl, ppm guess for 1H)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `kfigure()`, `kxlabel()`, `time()`, `kylabel()`, `xlim()`, `ylim()`, `xix_rep_time_ensemble_r()`, `T2n()`, `klegend()`, `savefig()`, `gaussleg()`, `r1n_dnp()`, `create()`, `basis()`, `state()`, `rep_time()`.
