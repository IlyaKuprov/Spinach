# examples/nmr_solids/cp_powder_static_nhh.m

https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/cp_powder_static_nhh.m

## Purpose

Static-powder 1H–15N cross-polarisation in the doubly rotating frame for one 15N surrounded by eight protons. The source describes the proton bath as scattered on a 2 Å-radius sphere around the nitrogen and gives an estimated runtime of minutes on a Tesla A100, much longer on CPU; these are source notes, not a recorded benchmark from this task.

## Spin model and reduced basis

The code sets the field parameter to 9.394, gives coordinates for all nine spins, and assigns scalar Zeeman values individually. The coordinate geometry defines the proton–nitrogen dipolar network. It requests anisotropic equilibrium at 298 K. The sphten-liouv basis uses IK-0 with inter_level 4; the source describes the retained space as including correlations through four spins. The source comment says a GPU is needed and the code enables greedy mode.

## Cross-polarisation observable

Both irradiation-power rows contain 5e4 over 100 points; the source gives no unit for this value. The RF operators are Ly on 1H and Lx on 15N, with excitation operators Lx on 1H and Ly on 15N. Detection is the 15N Lx state. The powder grid is rep_2ang_100pts_sph, and 100 steps of 1e-5 seconds span 1 ms. The wrapper calls cp_contact_hard and plots the real 15N response versus cumulative time. It does not specify a Hartmann–Hahn matching condition or expose the helper's internal contact dynamics. No measured spectrum is reported.
