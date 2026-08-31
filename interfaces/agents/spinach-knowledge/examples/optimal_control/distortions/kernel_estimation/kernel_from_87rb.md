# examples/optimal_control/distortions/kernel_estimation/kernel_from_87rb.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/distortions/kernel_estimation/kernel_from_87rb.m`
- Signature: `kernel_from_87rb()`
- Total lines: 130

## Purpose

Transmitter and probe distortion kernel of a 400 MHz Bruker spectrometer fitted with a 4 mm Phoenix MAS probe, estimated from an oscilloscope recording of an eight-block XiX wave- form on the 87Rb channel. The wall clock record is heterodyned into the rotating frame, the carrier frequency is refined from the residual phase drift accumulated within each XiX half-block, the remaining constant phase is removed, and the 

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 22-23: Kernel sampling parameters; implemented by `dt_kernel=1e-7; n_taps=40; t_kernel=dt_kernel*n_taps`.
- Lines 25-26: XiX waveform parameters; implemented by `n_blocks=8; t_block=10e-6`.
- Lines 28-29: Start of the XiX waveform in the oscilloscope record; implemented by `t_start=1.11e-5`.
- Lines 31-32: Nominal carrier frequency of the 87Rb channel; implemented by `carrier_freq=131.1468520e6`.
- Lines 34-35: Fractional window inside each half-block, clear of the transitions; implemented by `win_start=0.32; win_end=0.48`.
- Lines 37-38: Read the oscilloscope record; implemented by `load('87Rb_160W_XiX_10us_8repeats.mat','scope_trace','time_scope')`.
- Lines 41-42: Plot the raw scope readout; implemented by `kfigure(); scale_figure([2.0 1.5])`.
- Lines 48-49: Segment holding the waveform and its ringdown; implemented by `n_start=round(t_start/dt)`.
- Lines 52-53: Time axis of the segment, relative to the start of the waveform; implemented by `time_grid=dt*(0:(n_end-n_start))'`.
- Lines 55-56: Fitting windows inside each XiX half-block; implemented by `n_low=round((win_start+0.5*(0:(2*n_blocks-1))')*t_block/dt)+1`.
- Lines 59-60: Heterodyne the record at the nominal carrier frequency; implemented by `[real_bruk,imag_bruk]=heterodyne(dt,scope_trace,carrier_freq)`.
- Lines 63-64: Residual carrier offset from the phase slope in each half-block; implemented by `freq_drift=zeros(2*n_blocks,1)`.
- Lines 72-73: Refine the carrier frequency, skipping the leading transient block; implemented by `carrier_freq=carrier_freq-mean(freq_drift(2:end))/(2*pi)`.
- Lines 75-76: Heterodyne the record again at the refined carrier frequency; implemented by `[real_bruk,imag_bruk]=heterodyne(dt,scope_trace,carrier_freq)`.
- Lines 79-80: Constant phase of the nominally positive XiX half-blocks; implemented by `block_phase=zeros(n_blocks,1)`.
- Lines 86-87: Remove the constant phase and normalise the amplitude; implemented by `cplx_bruk=exp(-1i*angle(mean(exp(1i*block_phase))))*cplx_bruk`.
- Lines 90-91: Resample onto the kernel time grid; implemented by `real_bruk=resample(real(cplx_bruk),round(1/dt_kernel),round(1/dt))`.
- Lines 96-97: Ideal XiX waveform, silent after the last block; implemented by `xix_ideal=(-1).^floor(time_grid/(0.5*t_block))`.

### Control flow inferred from the code

- Line 65: `for` loop over `n=1:(2*n_blocks)`.
- Line 81: `for` loop over `n=1:n_blocks`.

### Key state/data transformations

- Lines 23: computes `dt_kernel` using `dt_kernel=1e-7; n_taps=40; t_kernel=dt_kernel*n_taps`.
- Lines 26: computes `n_blocks` using `n_blocks=8; t_block=10e-6`.
- Lines 29: computes `t_start` using `t_start=1.11e-5`.
- Lines 32: computes `carrier_freq` using `carrier_freq=131.1468520e6`.
- Lines 35: computes `win_start` using `win_start=0.32; win_end=0.48`.
- Lines 39: computes `dt` using `dt=time_scope(2)-time_scope(1)`.
- Lines 49: computes `n_start` using `n_start=round(t_start/dt)`.
- Lines 50: computes `n_end` using `n_end=n_start+round((n_blocks*t_block+t_kernel)/dt)`.
- Lines 53: computes `time_grid` using `time_grid=dt*(0:(n_end-n_start))'`.
- Lines 56: computes `n_low` using `n_low=round((win_start+0.5*(0:(2*n_blocks-1))')*t_block/dt)+1`.
- Lines 57: computes `n_up` using `n_up=round((win_end+0.5*(0:(2*n_blocks-1))')*t_block/dt)+1`.
- Lines 60: computes `[real_bruk,imag_bruk]` using `[real_bruk,imag_bruk]=heterodyne(dt,scope_trace,carrier_freq)`.
- Lines 61: computes `cplx_bruk` using `cplx_bruk=real_bruk(n_start:n_end)+1i*imag_bruk(n_start:n_end)`.
- Lines 64: computes `freq_drift` using `freq_drift=zeros(2*n_blocks,1)`.
- Lines 66: computes `fit_win` using `fit_win=n_low(n):n_up(n)`.
- Lines 67-68: computes `phase_fit` using `phase_fit=polyfit(time_grid(fit_win), unwrap(angle(cplx_bruk(fit_win))),1)`.
- Lines 69: computes `freq_drift(n)` using `freq_drift(n)=phase_fit(1)`.
- Lines 80: computes `block_phase` using `block_phase=zeros(n_blocks,1)`.

## Implementation structure

- Transmitter and probe distortion kernel of a 400 MHz Bruker
- spectrometer fitted with a 4 mm Phoenix MAS probe, estimated
- from an oscilloscope recording of an eight-block XiX wave-
- form on the 87Rb channel.
- The wall clock record is heterodyned into the rotating frame,
- the carrier frequency is refined from the residual phase drift
- accumulated within each XiX half-block, the remaining constant
- phase is removed, and the FIR kernel is then obtained by line-
- ar least squares from the ideal and the measured waveforms.
- The heterodyne is an analytic signal demodulation, which is
- zero-phase: the rotating frame signal is not delayed with re-
- spect to the oscilloscope record.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `load()`, `time_scope()`, `kfigure()`, `scale_figure()`, `subplot()`, `kxlabel()`, `kylabel()`, `ktitle()`, `heterodyne()`, `real_bruk()`, `imag_bruk()`, `n_low()`, `n_up()`, `polyfit()`, `time_grid()`, `unwrap()`.
