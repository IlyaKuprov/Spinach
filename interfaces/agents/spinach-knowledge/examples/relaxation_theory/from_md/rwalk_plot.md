# examples/relaxation_theory/from_md/rwalk_plot.m

- MATLAB implementation: [examples/relaxation_theory/from_md/rwalk_plot.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/from_md/rwalk_plot.m)

- Signature: `rwalk_plot()`
- Run from MATLAB with no arguments. The function creates a 3-D plot in the current graphics context and returns no MATLAB value.

## Purpose

Visualises one sampled rotational random walk on the unit sphere. It is a trajectory illustration, not a relaxation-rate calculation or a statistical validation of the walk.

## Sampling and plotting

The example sets `tau_c=6e-9`, takes `tc_steps=500` steps per correlation-time interval, and therefore uses `dt=tau_c/tc_steps`; it requests `tot_nsteps=5000` samples from `rwalk(tot_nsteps,tau_c,dt)`. It first calls `rng('shuffle')`, so successive runs are intentionally not reproducible from a fixed seed. The source does not state units for `tau_c` or `dt`; preserve the numeric values and the unit convention used by `rwalk`.

Each returned Euler-angle triple is converted by `euler2dcm` and applied to `[0;0;1]`, giving one point on the unit sphere. The points are connected with `plot3`; all three axes are limited to `[-1,1]`, with square aspect, grid, and box enabled. The example does not return the angles or trajectory, expose sampling controls as function arguments, or label the plot with time or step indices.
