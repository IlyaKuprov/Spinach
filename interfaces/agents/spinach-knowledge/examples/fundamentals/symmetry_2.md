# examples/fundamentals/symmetry_2.m

- MATLAB implementation: [examples/fundamentals/symmetry_2.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/symmetry_2.m)

- Signature: `symmetry_2()`
- Source: [examples/fundamentals/symmetry_2.m](../../../../../examples/fundamentals/symmetry_2.m)

## Purpose

Simulates a liquid-state pulse-acquire proton NMR FID for a highly symmetric 13-proton spin system, then apodises, Fourier-transforms, and plots its real spectrum. The source credits the spin system to Andres Castillo.

## Spin model and basis

The field is set to 9.4 T. The 13 proton Zeeman scalar values, in spin order, are 0.89, 0.89, 0.89, 0.895, 0.895, 0.895, 1.16, 1.16, 1.16, 1.2, 1.39, 1.71, and 3.85. The scalar-coupling matrix is divided by 2 in the source; its nonzero values before that division are 6.7200, 6.6400, 6.0800, 14.0850, 5.1000, 8.3000, 8.3050, and 6.4000.

The spherical-tensor Liouville basis uses the IK-2 approximation and three S3 groups on spins [1 2 3], [4 5 6], and [7 8 9], with scalar-coupling connectivity, proximity level 1, and projection +1. The source describes this as the fully symmetric irreducible representation of S3(x)S3(x)S3.

## Acquisition and processing

The initial state uses `state(spin_system,'L+','1H')`; no decoupled spins are specified. `liquid(...,@acquire,...,'nmr')` generates a 2048-point FID with sweep 2000 Hz and offset 800 Hz. The FID is apodised with `{'exp',6}`, zero-filled to 8196 points, Fourier-transformed with `fftshift(fft(...))`, and plotted as its real part with axis units set to ppm and the axis inverted. The receiver uses the same operator description with `coil_state` instead.

## Scope

The file specifies a simulation and plotting workflow; this entry does not assert a particular spectrum, peak position, or successful numerical run.
