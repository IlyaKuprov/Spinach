# examples/optimal_control/distortions/kernel_estimation/kernel_from_antenna.m

- MATLAB implementation: [examples/optimal_control/distortions/kernel_estimation/kernel_from_antenna.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/distortions/kernel_estimation/kernel_from_antenna.m)

Source: [examples/optimal_control/distortions/kernel_estimation/kernel_from_antenna.m](../../../../../../../examples/optimal_control/distortions/kernel_estimation/kernel_from_antenna.m)

- Signature: `kernel_from_antenna()`

## Purpose

Estimates a causal HiPER instrument-filter kernel by fitting a synthetic ideal XiX reference to the measured in-phase and quadrature voltage recorded by a nearby antenna. The method uses a previously saved acquisition; it does not acquire or validate a hardware response.

## Inputs and reference waveform

The script loads `time_ns`, `real_part`, and `imag_part` from `xix_on_resonance.mat`. These are the imported antenna record, plotted over 60–500 ns as real and imaginary channels in mV. The acquisition file and its measurement conditions are not specified by the script.

The reference is synthetic: each XiX block concatenates 36 samples at +1 and 36 at −1, and ten such blocks are placed between zero pads. A two-sample shift is applied; the real reference has 1152 samples in total and its quadrature component is zero. The source labels this construction as ten 36 ns periods. It is plotted on the supplied `time_ns` axis; the script does not state a separate sample interval.

## Kernel fit

The real and imaginary reference components form the complex input `x`; the two imported antenna channels form complex output `y`. The script calls `kernelest(x,y,32,'tikh','causal',10)`: a 32-point causal kernel is fitted using Tikhonov regularisation with parameter 10. The resulting `h` is saved to `hiper_kernel_antenna.mat`.

## Outputs and limits

The figure compares the measured antenna record with the ideal waveform, then plots the kernel's real and imaginary components against the helper-generated time axis. A second plot shows the magnitude of a 128-point shifted FFT, with its frequency axis constructed by `fft_freq_axis(32,0.5,96)` and labelled in GHz. These are plots of the fit and its derived frequency response, not an independent hardware-validation test; the script supplies no uncertainty estimate or acquisition provenance.

Contributors named in the source: Rob Hunter, Hassane el-Mkami, Graham Smith, Yujie Zhao, Shebha Anandhi Jegadeesan, Guinevere Mathies, and Ilya Kuprov.
