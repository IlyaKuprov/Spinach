# examples/optimal_control/distortions/rlc_response_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/distortions/rlc_response_1.m`
- Signature: `rlc_response_1(interp_type)`
- Total lines: 80

## Purpose

An illustration of the effect of the resonator response function on a typical composite pulse in NMR spectroscopy. The argument may be set to 'previous' (default, corresponds to piecewise con- stant input waveform) or any of the options ('linear', 'cubic', etc.) supported by interp1() function. Calculation time: seconds.

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.

## Numerical / algorithmic content

## Implementation structure

- An illustration of the effect of the resonator response function
- on a typical composite pulse in NMR spectroscopy. The argument
- may be set to 'previous' (default, corresponds to piecewise con-
- stant input waveform) or any of the options ('linear', 'cubic',
- etc.) supported by interp1() function.
- Calculation time: seconds.
- Default to piecewise-constant
- Decide the time grid (4 x Nyquist)
- Random amplitude component
- Random phase component
- Put the pulse together
- Heterodyne out the carrier frequency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `spin()`, `exist()`, `amp_part()`, `phi_part()`, `kfigure()`, `scale_figure()`, `subplot()`, `kxlabel()`, `kylabel()`, `ktitle()`, `ylim()`, `klegend()`, `lsim()`, `lowpass()`.
