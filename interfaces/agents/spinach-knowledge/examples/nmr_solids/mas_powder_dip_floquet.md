# examples/nmr_solids/mas_powder_dip_floquet.m

Source: [examples/nmr_solids/mas_powder_dip_floquet.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/mas_powder_dip_floquet.m)

## Model

This is a simulated MAS powder pulse-acquire example for two protons. The input sets a 14.1 T field, isotropic Zeeman scalar values of 5.0 and -2.0, and positions [0, 0, 0] and [0, 3.9, 0.1]. The coordinates provide the geometry for the dipolar interaction when the Spinach system is built; the source does not set a separate dipolar-coupling constant. The coordinate and Zeeman scalar units are not stated in this file. No quadrupolar spin, NQI value, or explicit RF pulse or decoupling field is specified.

The calculation uses the spherical-tensor Liouville-space basis with no approximation and projection +1. It sets a rotor axis of [1, 1, 1], a rotor rate of 1000 Hz, maximum rank 17, the orientation grid `leb_2ang_rank_17`, and rotor-frame reference. Its Floquet acquisition is the call `floquet(spin_system,@acquire,parameters,'nmr')`; this is the algorithmic distinction from the companion `singlerot` Fokker-Planck and `gridfree` examples.

## Signal and display

The initial state and receiver coil are both the proton `L+` state. `@acquire` produces the FID using 512 points, a sweep setting of 2e4, offset 0, and zero filling to 4096. The display axis is configured in ppm with inversion enabled. The code applies exponential apodisation parameter 6, computes `fftshift(fft(fid,parameters.zerofill))`, and plots its real part with `plot_1d`. This plotted spectrum is a computed result from the configured model, not an experimentally measured spectrum.