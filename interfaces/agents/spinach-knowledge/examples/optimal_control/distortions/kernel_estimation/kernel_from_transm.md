# examples/optimal_control/distortions/kernel_estimation/kernel_from_transm.m

- Signature: `kernel_from_transm()`

## Purpose

Extracts the HiPER instrument response from a transmission measurement. An AWG linear chirp spans 93.3–94.7 GHz over 1400 ns; the recorded power is square-rooted to obtain amplitude, and the two spikes mark the chirp start and finish. Since the measurement has no phase information, the script constructs a minimum-phase kernel from the measured amplitude spectrum using the real-cepstrum (Kramers–Kronig) construction; a zero-phase alternative would be symmetric in time and non-causal. Rob Hunter, Hassane el-Mkami, Graham Smith, Yujie Zhao, Shebha Anandhi Jegadeesan, Guinevere Mathies, Ilya Kuprov.

## Method

The script clips negative measured powers, normalises the amplitude spectrum, applies a sin⁴ taper while leaving 93.4–94.6 GHz intact, and resamples the spectrum. It forms a real cepstrum with a −60 dB floor, folds it causally to obtain the minimum-phase response, then resamples the kernel on a 0.5 ns grid and saves 32 points as `hiper_kernel_trans.mat`.
