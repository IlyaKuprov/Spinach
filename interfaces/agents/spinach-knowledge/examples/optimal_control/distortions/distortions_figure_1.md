# examples/optimal_control/distortions/distortions_figure_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/distortions/distortions_figure_1.m`
- Signature: `distortions_figure_1()`
- Total lines: 73

## Purpose

Figure 1 from the paper by Rasulov and Kuprov:

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.

## Numerical / algorithmic content

## Implementation structure

- Figure 1 from the paper by Rasulov and Kuprov:
- Pulse sequence and its discretisation
- Convert waveform to mT
- Apply a cascade of two single-pole filters
- Original vs second-order low-pass filter
- Apply a cascade of three single-zero filters
- Original vs third-order high-pass filter
- Apply an RLC filter
- Original vs RLC filter

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `fapt2sfo()`, `spin()`, `spf()`, `kfigure()`, `scale_figure()`, `subplot()`, `wave()`, `wave_spf()`, `klegend()`, `kxlabel()`, `kylabel()`, `szf()`, `wave_szf()`.
