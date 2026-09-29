# examples/optimal_control/distortions/kernel_estimation/kernel_from_transm.m

- MATLAB implementation: [examples/optimal_control/distortions/kernel_estimation/kernel_from_transm.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/distortions/kernel_estimation/kernel_from_transm.m)

Source: [examples/optimal_control/distortions/kernel_estimation/kernel_from_transm.m](../../../../../../../examples/optimal_control/distortions/kernel_estimation/kernel_from_transm.m)

- Signature: `kernel_from_transm()`

## Purpose

Constructs a 32-sample minimum-phase estimate of the HiPER transmission-response kernel from imported power-versus-frequency data. Unlike the antenna route, this acquisition contains amplitude information only; the script supplies a minimum-phase reconstruction rather than measured phase.

## Input and spectral preparation

The script loads `freq_axis_ghz` and `power_at_eik` from `power_at_eik.mat`. The source describes the imported profile as the power response to an AWG linear chirp from 93.3 to 94.7 GHz over 1400 ns, with spikes marking the chirp start and finish. It is not a simulated response. Negative power entries are clipped to zero, then square-rooted to amplitude and normalised by the maximum.

A sin⁴ taper fades the spectrum outside the open interval 93.4–94.6 GHz while leaving samples inside that band intact. The weighted spectrum is downsampled by 1000 (the source calls `resample` with factors 1 and 1000), and the plot compares the original amplitude, taper, and filtered spectrum. The script then zero-fills the amplitude spectrum on both sides and shifts the 94.0 GHz centre to zero frequency.

## Minimum-phase kernel construction

The real cepstrum is computed from the log amplitude after flooring values at 1e−3 (−60 dB). A causal cepstral window doubles the positive-time terms, leaves the endpoints unchanged, and zeros the negative-time terms; transforming back yields the minimum-phase response. The source uses the measured magnitude profile and this phase assumption, not a separately recorded phase trace.

The inverse-transform time axis is derived from the resampled frequency spacing. Samples through 16 ns are retained and scaled by 0.5 ns divided by the original time step; spline interpolation places the result on a 0.5 ns grid from 0 to 15.5 ns. The 32-point complex kernel is saved as `hiper_kernel_trans.mat`.

## Outputs and limits

The figure shows the input transmission amplitude with the taper and filtered curve, followed by the real and imaginary kernel components versus time. The latter is titled as the kernel at a 94.0 GHz offset. No spin-control pulse, field setting, independently measured phase, or hardware-validation run is part of this script.

Contributors named in the source: Rob Hunter, Hassane el-Mkami, Graham Smith, Yujie Zhao, Shebha Anandhi Jegadeesan, Guinevere Mathies, and Ilya Kuprov.
