# examples/optimal_control/distortions/kernel_estimation/kernel_from_transm.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/distortions/kernel_estimation/kernel_from_transm.m`
- Signature: `kernel_from_transm()`
- Total lines: 57

## Purpose

Response function extraction from the transmission profile of the HiPER instrument. We are sending a linear chirp via the AWG from 93.3 GHz to 94.7 GHz over 1400 ns and record- ing the power response (so a square root must be taken to obtain the amplitude). The two spikes show when the chirp starts and finishes. Rob Hunter, Hassane el-Mkami, Graham Smith, Yujie Zhao, Shebha Anandhi Jegadeesan, Guinevere Mathies, Ilya

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-16: Load power data and convert to amplitude; implemented by `load('power_at_eik.mat','freq_axis_ghz', 'power_at_eik')`.
- Lines 20-21: Plot as received; implemented by `kfigure(); scale_figure([0.75 1.00])`.
- Lines 27-29: Get apodisation weights; implemented by `leave_intact=(freq_axis_ghz>93.4)& (freq_axis_ghz<94.6)`.
- Lines 34-35: Apply weights and resample; implemented by `amp=resample(amp.*weights',1,1000)`.
- Lines 43-44: Get the convolution kernel with a 0.5 ns time step; implemented by `h=ifft(ifftshift([zeros(size(amp)); zeros(size(amp)); amp`.

### Key state/data transformations

- Lines 17: computes `power_at_eik(power_at_eik<0)` using `power_at_eik(power_at_eik<0)=0`.
- Lines 18: computes `amp` using `amp=sqrt(power_at_eik); amp=amp/max(amp)`.
- Lines 28-29: computes `leave_intact` using `leave_intact=(freq_axis_ghz>93.4)& (freq_axis_ghz<94.6)`.
- Lines 30: computes `np_intact` using `np_intact=nnz(leave_intact)`.
- Lines 31: computes `fade_in` using `fade_in=sin(linspace(0,pi/2,nnz(~leave_intact)/2)).^4`.
- Lines 32: computes `weights` using `weights=[fade_in ones(1,np_intact) fliplr(fade_in)]`.
- Lines 36-37: computes `freq_axis_ghz` using `freq_axis_ghz=linspace(freq_axis_ghz(1), freq_axis_ghz(end),numel(amp))'`.
- Lines 44: computes `h` using `h=ifft(ifftshift([zeros(size(amp)); zeros(size(amp)); amp`.
- Lines 46: computes `df` using `df=freq_axis_ghz(2)-freq_axis_ghz(1); zf=2*numel(amp)`.
- Lines 47: computes `[~,t]` using `[~,t]=ifft_time_axis(numel(amp),df,zf)`.

## Implementation structure

- Response function extraction from the transmission profile
- of the HiPER instrument. We are sending a linear chirp via
- the AWG from 93.3 GHz to 94.7 GHz over 1400 ns and record-
- ing the power response (so a square root must be taken to
- obtain the amplitude). The two spikes show when the chirp
- starts and finishes.
- Rob Hunter, Hassane el-Mkami, Graham Smith,
- Yujie Zhao, Shebha Anandhi Jegadeesan,
- Guinevere Mathies, Ilya Kuprov
- Load power data and convert to amplitude
- Plot as received
- Get apodisation weights

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `load()`, `power_at_eik()`, `kfigure()`, `scale_figure()`, `subplot()`, `kxlabel()`, `kylabel()`, `ktitle()`, `nnz()`, `fliplr()`, `resample()`, `freq_axis_ghz()`, `weights()`, `klegend()`, `ifftshift()`, `ifft_time_axis()`.
