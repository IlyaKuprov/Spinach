# examples/imaging/slice_select_1d_shaped.m

- MATLAB implementation: [examples/imaging/slice_select_1d_shaped.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/imaging/slice_select_1d_shaped.m)

## Purpose

This zero-argument example models one-dimensional slice selection with a Gaussian-shaped RF pulse while diffusion and flow are included in the imaging parameters. The source estimates seconds of calculation time.

## Spin, relaxation, and spatial model

The spin system is one 1H at 5.9 T with zero chemical shift. It uses diagonal T1/T2 relaxation, zero equilibrium, and rate values r1=30 and r2=70; the example does not attach units to these rates. The basis uses sphten-liouv without an approximation.

The sample extent is 0.30 on a 500-point grid, and the signal uses 128 points. The example sets both slice-selection and readout gradient amplitudes to 30e-3; the sequence callback documentation specifies T/m for the slice-selection gradient. The source does not label the sample-extent unit. The initial-state and receive-coil spatial profiles are uniform, with Lz initial state and L+ detection for 1H.

## Shaped pulse, transport, and readout

The RF pulse is represented by 50 steps over a total duration of 0.5e-4 s. Its frequency is +100e3 Hz, phase is pi/2, and per-step amplitudes are 2*pi*20000 multiplied by the Gaussian profile returned by pulse_shape. The sequence callback documentation specifies RF amplitude in rad/s and pulse duration in seconds; each step receives duration 0.5e-4/50 s. Maximum rank is 3. The flow-field vector is set to 1e-2 at every point and the diffusion parameter is 5e-6; the example does not state units for either value.

The example calls imaging with slice_select_1d in the imaging context, applies square-sine apodisation, computes a real shifted Fourier transform, and plots the resulting one-dimensional profile. It supplies transport parameters but does not state a numerical simulation outcome.

## Source

[Example source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/imaging/slice_select_1d_shaped.m) · [Imaging callback](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/imaging/slice_select_1d.m)

The source credits Ahmed Allami and Ilya Kuprov.
