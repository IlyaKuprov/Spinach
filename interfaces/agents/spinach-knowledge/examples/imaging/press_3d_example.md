# examples/imaging/press_3d_example.m

- MATLAB implementation: [examples/imaging/press_3d_example.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/imaging/press_3d_example.m)

## Purpose

This zero-argument example computes and plots a three-dimensional PRESS excitation profile using a tilted gradient system. Its comments say that changing pulse frequencies moves the hot spot through the sample. The source estimates hours of calculation time and notes that a Tesla V100 GPU makes it faster.

## Spin and spatial model

The 3.0 T model contains two 1H spins with Zeeman scalar values 1.5 and 3.7 and scalar coupling value 10; the example does not state units for these values. It uses the sphten-liouv formalism without a basis approximation, disables path tracing, and enables the polyadic option (the adjacent GPU option is commented out).

The uniform-density phantom has dimensions [0.30 0.25 0.27] and a [108 90 111] point spatial grid; the image grid is [101 103 105]. Three slice-selection gradients are set to 25e-3, and the spatial-selection gradient is 5e-3 T/m for 5e-4 s. The three gradient angles are [pi/3, pi/4, pi/5]. The plotted field-of-view axes are labelled in metres, while the example does not separately annotate the units of its dimensions parameter.

## Excitation-profile calculation

The RF phases are all pi/2; frequencies are {-30e3, +30e3, -50e3}, amplitudes are 2*pi*5000, durations are {0.5e-4, 1.0e-4, 1.0e-4} s, and maximum ranks are {2, 2, 2}. The imaging callback documentation gives Hz for RF frequency, rad/s for RF amplitude, and seconds for duration. A uniform phantom, total Lz initial state, and uniform receive-coil phantom are provided.

Unlike the 1D and 2D examples, this file calls imaging with press_voxel_3d only. It plots the resulting three-dimensional excitation profile with volplot; it does not run a PRESS spectral acquisition, Fourier transform, or spectrum plot. The example does not report a numerical excitation-profile result.

## Source

[Example source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/imaging/press_3d_example.m) · [Voxel diagnostic](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/imaging/press_voxel_3d.m)

The source credits Ahmed Allami and Ilya Kuprov.
