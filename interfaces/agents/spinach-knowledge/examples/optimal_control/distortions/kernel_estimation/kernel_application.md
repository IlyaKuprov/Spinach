# examples/optimal_control/distortions/kernel_estimation/kernel_application.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/distortions/kernel_estimation/kernel_application.m`
- Signature: `kernel_application()`
- Total lines: 44

## Purpose

HiPER instrument filter function kernel application to a complicated shaped pulse and a comparison with expe- rimental measurement at the instrument. Rob Hunter, Hassane el-Mkami, Graham Smith, Yujie Zhao, Shebha Anandhi Jegadeesan, Guinevere Mathies, Ilya Kuprov

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-13: Read the input pulse; implemented by `load('shaped_pulse_inp.mat','time_ns', 'real_part','imag_part')`.
- Lines 15-16: Plot the input pulse; implemented by `kfigure(); scale_figure([1.0 1.5])`.
- Lines 21-22: Load appropriate HiPER kernel; implemented by `load('hiper_kernel_antenna.mat','h')`.
- Lines 24-25: Compute and plot the convolution; implemented by `x=real_part+1i*imag_part`.
- Lines 31-33: Read the measured pulse; implemented by `load('shaped_pulse_out.mat','time_ns', 'real_part','imag_part')`.
- Lines 39-40: Plot the measured pulse; implemented by `subplot(3,1,3); plot(time_ns,[real_part imag_part])`.

### Key state/data transformations

- Lines 25: computes `x` using `x=real_part+1i*imag_part`.
- Lines 26: computes `y` using `y=conv(x,h); y=y(1:numel(x))`.
- Lines 34: computes `real_part` using `real_part=real_part(174:934)`.
- Lines 35: computes `imag_part` using `imag_part=imag_part(174:934)`.
- Lines 36: computes `time_ns` using `time_ns=time_ns(174:934)`.

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
