# examples/optimal_control/distortions/distortions_figure_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/distortions/distortions_figure_1.m`
- Signature: `distortions_figure_1()`
- Total lines: 73

## Purpose

Figure 1 from the paper by Rasulov and Kuprov:

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Pulse sequence and its discretisation; implemented by `fapt={[0 pi/10e-6 0 0.0e-6 5.0e-6]`.
- Lines 19-20: Convert waveform to mT; implemented by `wave=1e3*wave/spin('1H')`.
- Lines 22-23: Apply a cascade of two single-pole filters; implemented by `wave_spf=spf(wave, 0.9*exp(-1i*0.05))`.
- Lines 26-27: Original vs second-order low-pass filter; implemented by `kfigure(); scale_figure([1 1.6])`.
- Lines 39-40: Apply a cascade of three single-zero filters; implemented by `wave_szf=szf(wave, 0.1*exp(-1i*0.05))`.
- Lines 44-45: Original vs third-order high-pass filter; implemented by `subplot(3,1,2); hold on; box on`.
- Lines 56-57: Apply an RLC filter; implemented by `wave_spf=spf(wave, 0.9*exp(-1i*0.5))`.
- Lines 60-61: Original vs RLC filter; implemented by `subplot(3,1,3); hold on; box on`.

### Key state/data transformations

- Lines 11: computes `fapt` using `fapt={[0 pi/10e-6 0 0.0e-6 5.0e-6]`.
- Lines 16: computes `time_grid` using `time_grid=linspace(-5e-6,40e-6,1000)`.
- Lines 17: computes `wave` using `wave=fapt2sfo(fapt,time_grid)`.
- Lines 23: computes `wave_spf` using `wave_spf=spf(wave, 0.9*exp(-1i*0.05))`.
- Lines 40: computes `wave_szf` using `wave_szf=szf(wave, 0.1*exp(-1i*0.05))`.

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
