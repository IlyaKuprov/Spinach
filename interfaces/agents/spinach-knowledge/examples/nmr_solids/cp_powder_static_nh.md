# examples/nmr_solids/cp_powder_static_nh.m

https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/cp_powder_static_nh.m

## Purpose

A two-spin 1H–15N cross-polarisation calculation in the doubly rotating frame, averaged over a static powder. The source estimates seconds of calculation time; this is not a reported experimental or validation result.

## Spin model and experiment

The source sets sys.magnet to 9.394, places 1H and 15N at [0, 0, 0] and [0, 0, 1.05], and sets both isotropic Zeeman terms to zero. The coordinates provide the pair geometry for the dipolar interaction; the wrapper assigns no additional interaction tensor. It uses the full sphten-liouv basis with no basis approximation and requests anisotropic equilibrium at 298 K. The wrapper passes two irradiation-power rows, each 5e4 at 100 points, with Ly on 1H and Lx on 15N; the power units are not stated. The excitation operators are Lx on 1H and Ly on 15N, and the detected 15N state is Lx.

## Powder signal

The static powder calculation uses rep_2ang_6400pts_sph and 100 time steps of 1e-5 seconds (1 ms in total) with cp_contact_hard. The wrapper supplies those irradiation and time-grid parameters but does not state a Hartmann–Hahn matching criterion or define the helper's contact dynamics; neither should be inferred from the power values alone. The example plots the real 15N FID against cumulative time and labels it as the 15N S_X expectation value. No measured spectrum is reported.
