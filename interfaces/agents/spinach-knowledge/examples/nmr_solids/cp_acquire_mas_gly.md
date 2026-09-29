# examples/nmr_solids/cp_acquire_mas_gly.m

Source: [examples/nmr_solids/cp_acquire_mas_gly.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/cp_acquire_mas_gly.m)

## Purpose and spin model

This example simulates 1H-to-13C cross-polarisation (CP), followed by acquisition under MAS in alpha-glycine powder. Spin-system properties are read from a Gaussian log parsed with `gparse` and passed to `g2spinach`; the source labels these properties as PCM DFT. It then labels the spectrometer 400 MHz, sets `sys.magnet=9.4`, and assigns alpha-glycine isotropic shifts: CO 176.4, CA 43.6, the two H_CA sites 2.6 and 3.8, and three H_N sites 8.0. The source does not label units for these shift values. The spin temperature parameter is 298, with no unit stated in this file.

The basis is `sphten-liouv` with the `IK-0` approximation and `inter_level=3`, retaining correlations up to and including three-spin correlations as stated in the header. The source enables the greedy option, disables `pt`, and sets an interaction cutoff expression `2*pi*200`; its comment describes the cutoff as neglecting interactions below 200 Hz. No claim about a measured spectrum is made by these settings.

## CP, MAS, and readout

The example sets `parameters.spins={'1H','13C'}`, `rate=10000`, axis `[sqrt(2/3) 0 sqrt(1/3)]`, maximum rank 5, and powder grid `rep_2ang_100pts_sph`. The source does not attach a unit to the rotor-rate value. It specifies offsets `[2e3 10e3]` Hz, high-power value `83e3` Hz, CP powers `[60e3 50e3]` Hz, CP duration `50e-5` seconds, and acquisition sweep `5e4` Hz. Acquisition uses 512 points and zero-fills to 4096. The source requests the `iso_eq` initial condition and uses the 13C `L+` state as the receive coil.

The simulation calls `singlerot` with `@cp_acquire_soft` in NMR mode, applies exponential apodisation with parameter 6, Fourier transforms the FID, and plots its real part for 13C using the second configured offset. Thus the plotted result is a simulated carbon spectrum after CP and MAS, not an experimental input spectrum. The source reports minutes on a Tesla A100 and much longer on CPU; it gives no peak values or comparison to measured data.
