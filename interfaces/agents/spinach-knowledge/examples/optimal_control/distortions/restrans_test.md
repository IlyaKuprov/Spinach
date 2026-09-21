# examples/optimal_control/distortions/restrans_test.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/distortions/restrans_test.m`
- Signature: `restrans_test()`
- Total lines: 38

## Purpose

Resonator transform test. Sends a square pulse into a simple resonator model and plots the time-domain response.

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.

## Numerical / algorithmic content

## Implementation structure

- Resonator transform test. Sends a square pulse into a simple
- resonator model and plots the time-domain response.
- Get a pulse shape
- Get a figure going

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `pulse_in()`, `kfigure()`, `scale_figure()`, `subplot()`, `restrans()`, `ktitle()`.
