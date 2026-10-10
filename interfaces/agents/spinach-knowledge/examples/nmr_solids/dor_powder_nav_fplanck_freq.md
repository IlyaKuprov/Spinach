# examples/nmr_solids/dor_powder_nav_fplanck_freq.m

https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/dor_powder_nav_fplanck_freq.m

## Purpose and model

A frequency-domain double-angle-spinning calculation for the 14N nucleus of N-acetylvaline. The source says the 1D Fokker–Planck treatment includes the second-order quadrupolar shift and third-order lineshape. It warns that the outer and inner spinning frequencies are set artificially high to shorten this example, and that slower rates or larger NQIs need larger ranks and spherical grids. The stated calculation time is minutes.

The one-spin model uses sys.magnet=14.1 and constructs its quadrupolar coupling through eeqq2nqi(3.21e6, 0.27, 1, [0, 0, 0]); the source does not attach units to these arguments. Relaxation is diagonal damping with zero equilibrium and damp_rate=2e3 (unit not stated). The full sphten-liouv basis has no approximation.

## DOR and observable

The wrapper sets rate_outer=1e6 and rate_inner=5e6, with ranks 7 and 4. Its axis comments identify 54.74° for the outer rotor and 30.56° for the inner rotor. It uses rep_2ang_100pts_oct, a sweep parameter of [-50000, 50000], 1024 points, 1024-point zero filling, and kHz as the axis unit. Initial state and coil are both 14N L+, with the 14N rotating frame set to 3. The calculation calls doublerot with slowpass in the lab frame and plots the real frequency-domain spectrum. This wrapper does not define the slowpass sequence internals; its source parameters are not experimental measurements or a validation result.
