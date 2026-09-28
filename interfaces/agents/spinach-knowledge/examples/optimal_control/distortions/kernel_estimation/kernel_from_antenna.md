# examples/optimal_control/distortions/kernel_estimation/kernel_from_antenna.m

- Signature: `kernel_from_antenna()`

## Purpose

Estimates the HiPER instrument filter-function kernel from in-phase and quadrature components recorded by an antenna near the sample. Rob Hunter, Hassane el-Mkami, Graham Smith, Yujie Zhao, Shebha Anandhi Jegadeesan, Guinevere Mathies, Ilya Kuprov.

## Method

Loads `xix_on_resonance.mat`, plots the measured signal, and constructs an ideal XiX waveform with ten 36 ns periods and a 2 ns shift. It estimates a 32-point causal kernel using `kernelest` with Tikhonov regularisation and parameter 10, then plots the convolution and the kernel's time- and frequency-domain responses. The kernel is saved as `hiper_kernel_antenna.mat`.
