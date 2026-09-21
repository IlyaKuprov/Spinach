# examples/optimal_control/distortions/distortions_figure_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/distortions/distortions_figure_2.m`
- Signature: `distortions_figure_2()`
- Total lines: 39

## Purpose

Figure 2 from the paper by Rasulov and Kuprov:

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.

## Numerical / algorithmic content

## Implementation structure

- Figure 2 from the paper by Rasulov and Kuprov:
- Get E1000B pulse from Spinach
- Apply amplifier compression
- Plot waveforms
- Plot saturation levels
- Annotate the plot

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `vg_pulse()`, `amp_tanh()`, `tanh_waveform()`, `amp_root()`, `root_waveform()`, `kfigure()`, `kxlabel()`, `kylabel()`, `klegend()`.
