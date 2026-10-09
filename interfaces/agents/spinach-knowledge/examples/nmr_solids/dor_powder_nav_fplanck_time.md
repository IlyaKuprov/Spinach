# examples/nmr_solids/dor_powder_nav_fplanck_time.m

- Signature: dor_powder_nav_fplanck_time()

## Purpose

The example models a one-dimensional double-angle-spinning (DOR) powder spectrum for the ¹⁴N nucleus of N-acetylvaline. Its source describes time-domain detection, second-order quadrupolar shifts, and third-order lineshape contributions. The source comment estimates seconds of calculation and explicitly says its rotor rates are set artificially high to shorten this example; it warns that slower spinning or larger NQIs need larger ranks and spherical grids. These are source comments, not timing or convergence results from a run in this task.

## Spin model and rotor treatment

The system contains only ¹⁴N and one quadrupolar coupling tensor, supplied through eeqq2nqi(3.21e6, 0.27, 1, [0 0 0]). The source does not annotate the unit of 3.21e6. It sets sys.magnet to 14.1 without a unit comment. The basis is sphten-liouv with no approximation. Damping relaxation is selected with a code value of 2e3, diagonal relaxation retention, and zero equilibrium state.

The DOR propagation uses the one-dimensional Fokker–Planck method on the rep_2ang_100pts_oct spherical grid. Code-set outer and inner rotor-rate values are 1e6 and 5e6; the source does not state their units. The corresponding ranks are 7 and 4. The axis vectors are explicitly commented as 54.74° (outer) and 30.56° (inner). This is a single-spin powder calculation, not a CP or HMQC sequence.

## Acquisition and reported observable

The script sets a sweep of 1e5, 256 acquired points, and zero-fills to 1024; the displayed axis unit is kHz. Both initial state and receiver coil are ¹⁴N L+, and the rotating-frame entry is ¹⁴N, frame 3. It passes these settings and the acquire callback to doublerot in the lab frame, Fourier-transforms the resulting FID, and plots the real spectrum. The script defines that calculated observable. The wrapper call does not add an RF pulse, Hartmann–Hahn condition, or additional experimental internals.

Source: https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/dor_powder_nav_fplanck_time.m
