# examples/nmr_solids/mqmas_nqi.m

- MATLAB implementation: [examples/nmr_solids/mqmas_nqi.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/mqmas_nqi.m)

Source: [examples/nmr_solids/mqmas_nqi.m](../../../../../examples/nmr_solids/mqmas_nqi.m)

- Signature: `mqmas_nqi()`

## Model

The spin system contains `87Rb` at 9.4 T, with the source specifying just a nuclear quadrupole interaction. It constructs the NQI with `eeqq2nqi(5e6,0.50,3/2,[0 0 0])`: a 5 MHz coupling input, asymmetry parameter 0.50, spin 3/2, and zero Euler-angle values. These are settings for the simulated example, not a reported measurement of a particular compound.

## Rotor-synchronous MQMAS sequence

The script sets a 62.5 kHz rotor rate about axis vector `[1 1 1]`, rank 7, and powder grid `rep_2ang_1600pts_sph`. It uses MQ order 3, zero offset from the transmitter position described as the isotropic chemical shift, and puts `87Rb` in rotor frame 2. Two RF amplitudes are specified as `2*pi*250e3` rad/s each (250 kHz in cycles per second), with pulse durations 2 microseconds and 1 microsecond. The initial state is `Lz` and the receiver is `L+`.

The experiment is run with `singlerot` and `mqmas` in the lab frame. The two acquisition dimensions use 128 points each and are zero-filled to 256 each. Both dimensions receive squared-cosine apodisation before a two-dimensional Fourier transform; the magnitude spectrum is plotted. These settings define a calculation, not an experimentally measured spectrum.
