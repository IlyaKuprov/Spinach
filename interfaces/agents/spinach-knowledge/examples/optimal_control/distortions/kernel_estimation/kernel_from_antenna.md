# examples/optimal_control/distortions/kernel_estimation/kernel_from_antenna.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/optimal_control/distortions/kernel_estimation/kernel_from_antenna.m`
- Signature: `kernel_from_antenna()`
- Total lines: 61

## Purpose

HiPER instrument filter function kernel estimation from the quadrature components recorded by an antenna placed close to the sample location. Rob Hunter, Hassane el-Mkami, Graham Smith, Yujie Zhao, Shebha Anandhi Jegadeesan, Guinevere Mathies, Ilya Kuprov

## Physical / mathematical content

- Optimal-control examples. These scripts formulate pulse design as a nonlinear optimisation problem over waveform samples or basis coefficients. The core mathematical objects are fidelities, gradients, Hessians or Hessian approximations, ensemble robustness objectives, and constrained search over RF amplitude/phase trajectories.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Read the experimental data; implemented by `load('xix_on_resonance.mat','time_ns','real_part','imag_part')`.
- Lines 14-15: Plot the experimental data; implemented by `kfigure(); scale_figure([1.0 1.5])`.
- Lines 20-21: Make ideal XiX waveform with ten 36 ns periods; implemented by `xix_block=[ones(36,1); -ones(36,1)]; ts=2`.
- Lines 26-27: Plot the ideal waveform; implemented by `subplot(3,1,2); plot(time_ns,[real_ideal imag_ideal])`.
- Lines 31-32: Extract the kernel; implemented by `x=real_ideal+1i*imag_ideal`.
- Lines 37-38: Compute and plot the convolution; implemented by `y=conv(x,h); y=y(1:numel(x))`.
- Lines 43-44: Plot the kernel; implemented by `kfigure(); scale_figure([0.75 1.5])`.
- Lines 51-52: Plot the frequency response; implemented by `h=abs(fftshift(fft(h,128)))`.

### Key state/data transformations

- Lines 21: computes `xix_block` using `xix_block=[ones(36,1); -ones(36,1)]; ts=2`.
- Lines 22: computes `real_ideal` using `real_ideal=[zeros(175-ts,1); xix_block; xix_block; xix_block; xix_block; xix_block`.
- Lines 24: computes `imag_ideal` using `imag_ideal=zeros(1152,1)`.
- Lines 32: computes `x` using `x=real_ideal+1i*imag_ideal`.
- Lines 33: computes `y` using `y=real_part+1i*imag_part`.
- Lines 34: computes `h` using `h=kernelest(x,y,32,'tikh','causal',10)`.
- Lines 45: computes `time_axis` using `time_axis=0.5*(0:31)'; subplot(2,1,1)`.
- Lines 53: computes `freq_axis` using `freq_axis=fft_freq_axis(32,0.5,128-32)`.

## Implementation structure

- HiPER instrument filter function kernel estimation from
- the quadrature components recorded by an antenna placed
- close to the sample location.
- Rob Hunter, Hassane el-Mkami, Graham Smith,
- Yujie Zhao, Shebha Anandhi Jegadeesan,
- Guinevere Mathies, Ilya Kuprov
- Read the experimental data
- Plot the experimental data
- Make ideal XiX waveform with ten 36 ns periods
- Plot the ideal waveform
- Extract the kernel
- Compute and plot the convolution

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `load()`, `kfigure()`, `scale_figure()`, `subplot()`, `xlim()`, `kxlabel()`, `kylabel()`, `ktitle()`, `kernelest()`, `save()`, `fftshift()`, `fft_freq_axis()`.
