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
