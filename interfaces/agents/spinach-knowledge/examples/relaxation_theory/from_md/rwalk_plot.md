# examples/relaxation_theory/from_md/rwalk_plot.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/from_md/rwalk_plot.m`
- Signature: `rwalk_plot()`
- Total lines: 33

## Purpose

A plot of a typical random walk on a sphere.

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.

## Numerical / algorithmic content

## Implementation structure

- A plot of a typical random walk on a sphere.
- tau_c and number of steps per tau_c
- Total number of steps
- Euler angles of random walk
- Trajectory preallocation
- Trajectory
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `rng()`, `rwalk()`, `traj()`, `euler2dcm()`, `eulers()`, `plot3()`, `xlim()`, `ylim()`, `zlim()`.
