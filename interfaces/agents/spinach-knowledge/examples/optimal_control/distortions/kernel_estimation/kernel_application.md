# examples/optimal_control/distortions/kernel_estimation/kernel_application.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/distortions/kernel_estimation/kernel_application.m`
- Signature: `kernel_application()`
- Total lines: 44

## Purpose

HiPER instrument filter function kernel application to a complicated shaped pulse and a comparison with expe- rimental measurement at the instrument. Rob Hunter, Hassane el-Mkami, Graham Smith, Yujie Zhao, Shebha Anandhi Jegadeesan, Guinevere Mathies, Ilya Kuprov

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.

## Numerical / algorithmic content

## Implementation structure

- HiPER instrument filter function kernel application to
- a complicated shaped pulse and a comparison with expe-
- rimental measurement at the instrument.
- Rob Hunter, Hassane el-Mkami, Graham Smith,
- Yujie Zhao, Shebha Anandhi Jegadeesan,
- Guinevere Mathies, Ilya Kuprov
- Read the input pulse
- Plot the input pulse
- Load appropriate HiPER kernel
- Compute and plot the convolution
- Read the measured pulse
- Plot the measured pulse

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `load()`, `kfigure()`, `scale_figure()`, `subplot()`, `kxlabel()`, `kylabel()`, `ktitle()`, `real_part()`, `imag_part()`, `time_ns()`.
