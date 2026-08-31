# examples/optimal_control/pulse_analysis.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/pulse_analysis.m`
- Signature: `pulse_analysis()`
- Total lines: 23

## Purpose

An example of spectrogram analysis for a quadratic chirp pulse; adapted from Matlab example set. Calculation time: seconds.

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Quadratic chirp superposition; implemented by `fs=1000; t=0:1/fs:2-1/fs`.
- Lines 15-16: Do the plotting; implemented by `kfigure(); subplot(1,2,1); plot(t,y); kgrid`.

### Key state/data transformations

- Lines 11: computes `fs` using `fs=1000; t=0:1/fs:2-1/fs`.
- Lines 12-13: computes `y` using `y=chirp(t,100,1,200,'quadratic')+ chirp(t,200,1,100,'quadratic')`.

## Implementation structure

- An example of spectrogram analysis for a quadratic chirp
- pulse; adapted from Matlab example set.
- Calculation time: seconds.
- Quadratic chirp superposition
- Do the plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `chirp()`, `kfigure()`, `subplot()`, `kxlabel()`, `kylabel()`, `scale_figure()`, `spectrogram()`.
