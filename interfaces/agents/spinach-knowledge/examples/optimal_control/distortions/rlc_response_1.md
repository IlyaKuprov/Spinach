# examples/optimal_control/distortions/rlc_response_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/distortions/rlc_response_1.m`
- Signature: `rlc_response_1(interp_type)`
- Total lines: 80

## Purpose

An illustration of the effect of the resonator response function on a typical composite pulse in NMR spectroscopy. The argument may be set to 'previous' (default, corresponds to piecewise con- stant input waveform) or any of the options ('linear', 'cubic', etc.) supported by interp1() function. Calculation time: seconds.

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 18-19: Default to piecewise-constant; implemented by `if ~exist('interp_type','var')`.
- Lines 23-24: Decide the time grid (4 x Nyquist); implemented by `dt=pi/(4*omega); npts=tmax/dt+1`.
- Lines 27-28: Random amplitude component; implemented by `amp_part=rand(1,n_slices+1)`.
- Lines 33-34: Random phase component; implemented by `phi_part=2*pi*rand(1,n_slices+1)`.
- Lines 39-40: Put the pulse together; implemented by `inp_signal=amp_part.*cos(omega*time_grid+phi_part)`.
- Lines 42-43: Heterodyne out the carrier frequency; implemented by `inp_real=amp_part.*cos(phi_part)`.
- Lines 46-47: Plot input signal components; implemented by `kfigure(); scale_figure([1.75 1.2])`.
- Lines 58-59: Build the RLC bandpass response kernel; implemented by `sys=tf([1/(omega*Q) 0],[1/(omega^2) 1/(omega*Q) 1])`.
- Lines 61-62: Apply the RLC bandpass response kernel; implemented by `[out_signal,time_grid]=lsim(sys,inp_signal,time_grid)`.
- Lines 64-65: Heterodyne out the carrier frequency; implemented by `out_real=+2*lowpass(out_signal.*cos(omega*time_grid),0.1)`.
- Lines 68-69: Plot output signal components; implemented by `subplot(2,2,3); plot(1e6*time_grid,out_signal)`.

### Control flow inferred from the code

- Line 19: conditional branch on `~exist('interp_type','var')`.

### Key state/data transformations

- Lines 15: computes `omega` using `omega=14.09*spin('14N'); Q=80`.
- Lines 16: computes `n_slices` using `n_slices=20; tmax=50e-6`.
- Lines 20: computes `interp_type` using `interp_type='previous'`.
- Lines 24: computes `dt` using `dt=pi/(4*omega); npts=tmax/dt+1`.
- Lines 25: computes `time_grid` using `time_grid=linspace(0,tmax,npts)`.
- Lines 28: computes `amp_part` using `amp_part=rand(1,n_slices+1)`.
- Lines 29: computes `amp_part(1:2)` using `amp_part(1:2)=0; amp_part((end-2):end)=0`.
- Lines 34: computes `phi_part` using `phi_part=2*pi*rand(1,n_slices+1)`.
- Lines 35: computes `phi_part(1:2)` using `phi_part(1:2)=0; phi_part((end-2):end)=0`.
- Lines 40: computes `inp_signal` using `inp_signal=amp_part.*cos(omega*time_grid+phi_part)`.
- Lines 43: computes `inp_real` using `inp_real=amp_part.*cos(phi_part)`.
- Lines 44: computes `inp_imag` using `inp_imag=amp_part.*sin(phi_part)`.
- Lines 59: computes `sys` using `sys=tf([1/(omega*Q) 0],[1/(omega^2) 1/(omega*Q) 1])`.
- Lines 62: computes `[out_signal,time_grid]` using `[out_signal,time_grid]=lsim(sys,inp_signal,time_grid)`.
- Lines 65: computes `out_real` using `out_real=+2*lowpass(out_signal.*cos(omega*time_grid),0.1)`.
- Lines 66: computes `out_imag` using `out_imag=-2*lowpass(out_signal.*sin(omega*time_grid),0.1)`.

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
