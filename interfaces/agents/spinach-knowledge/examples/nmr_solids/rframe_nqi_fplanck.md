# examples/nmr_solids/rframe_nqi_fplanck.m

- Signature: `rframe_nqi_fplanck()`
- Source: [examples/nmr_solids/rframe_nqi_fplanck.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/rframe_nqi_fplanck.m)

## Purpose

Calculates a powder MAS spectrum with rotor-synchronised detection for a single quadrupolar 14N nucleus. The one-dimensional Fokker–Planck treatment uses a spherical grid and applies numerical second-order corrections to the rotating-frame transformation to account for the second-order quadrupolar shift and lineshape. The source estimates hours of calculation time.

## Spin system and rotor sampling

The field parameter is 14.1. The quadrupolar interaction is `eeqq2nqi(3.06e6, 0.40, 1, [0 0 0])`; its arguments have no units stated in the source. The full spherical-tensor Liouville basis is used. `singlerot` propagates the acquisition callback in the lab frame, with rate 50000, rotor axis `[1,1,1]`, maximum rank 85, and the `rep_2ang_6400pts_sph` grid. The example disables trajectory-level propagation and Krylov methods, and selects rotating-frame order 2. No gradient is configured.

## Acquisition and processing

Both the initial state and receiver are the 14N `L+` state. Acquisition uses sweep 50000, 256 points, zero-fill to 1024, and offset 18000; the frequency-axis units are explicitly set to Hz. The FID receives exponential apodisation with parameter 6, is Fourier transformed, and the real part is plotted.
