# examples/optimal_control/distortions/distortions_figure_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/distortions/distortions_figure_2.m`
- Signature: `distortions_figure_2()`
- Total lines: 39

## Purpose

Figure 2 from the paper by Rasulov and Kuprov:

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Get E1000B pulse from Spinach; implemented by `waveform=vg_pulse('E1000B',1000,0.001)`.
- Lines 14-15: Apply amplifier compression; implemented by `tanh_waveform=amp_tanh([waveform`.
- Lines 22-23: Plot waveforms; implemented by `kfigure(); hold on`.
- Lines 28-29: Plot saturation levels; implemented by `plot([0 1 NaN 1 0],[3e4 3e4 NaN -3e4 -3e4],'k--')`.
- Lines 31-32: Annotate the plot; implemented by `kxlabel('time, ms')`.

### Key state/data transformations

- Lines 11: computes `waveform` using `waveform=vg_pulse('E1000B',1000,0.001)`.
- Lines 12: computes `time_axis` using `time_axis=1e3*linspace(0,0.001,1000)`.
- Lines 15: computes `tanh_waveform` using `tanh_waveform=amp_tanh([waveform`.
- Lines 18: computes `root_waveform` using `root_waveform=amp_root([waveform`.
- Lines 34-35: computes `klegend({'input','output, tanh','output, $s` using `klegend({'input','output, tanh','output, $s=10$', 'saturation level'},'Location','Best')`.

## Implementation structure

- Figure 2 from the paper by Rasulov and Kuprov:
- Get E1000B pulse from Spinach
- Apply amplifier compression
- Plot waveforms
- Plot saturation levels
- Annotate the plot

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `vg_pulse()`, `amp_tanh()`, `tanh_waveform()`, `amp_root()`, `root_waveform()`, `kfigure()`, `kxlabel()`, `kylabel()`, `klegend()`.
