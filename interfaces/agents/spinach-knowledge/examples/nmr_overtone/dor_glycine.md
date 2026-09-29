# examples/nmr_overtone/dor_glycine.m

- MATLAB implementation: [examples/nmr_overtone/dor_glycine.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_overtone/dor_glycine.m)

This function sets up a simulated panoramic double-rotation (DOR) overtone spectrum for glycine's 14N nucleus. Its source comment says the simulation is described in Figure 1B of the paper at https://doi.org/10.1039/C5CP03266K, and attributes the glycine quadrupolar tensor data to O'Dell and Ratcliffe at https://doi.org/10.1016/j.cplett.2011.08.030. These are source attributions; the script does not itself report a fit, an experimental comparison, or validation against either paper. The comment estimates hours of calculation and says a short pulse with instrumentally inaccessible power is used to make the excitation pattern uniform.

The system has a single isotope, 14N, and a magnet setting of 14.1. The coupling input is exactly `eeqq2nqi(1.18e6,0.53,1,[0 0 0])`; the source does not annotate units for the magnet setting or these tensor arguments. Relaxation is `damp` with diagonal retention, zero equilibrium, and `damp_rate=100`. The basis is `sphten-liouv` with approximation `none`.

The DOR setup sets `theta=atan(sqrt(2))`, outer and inner rates 1425 and 6950, and rank 5 for each rotor. The axes are `[sin(theta) 0 cos(theta)]` and `[sqrt(20-2*sqrt(30)) 0 sqrt(15+2*sqrt(30))]`, commented as 54.74 and 30.56 degrees. The orientation grid is `rep_2ang_200pts_sph`. The spectrum settings are sweep `[-2e4 3e4]`, 1024 points, 1024-point zero-fill, and `axis_units='kHz'`. The initial state is 14N `Lz`; the receiver and Lx operator are the theta-weighted Lz/Lx combinations defined in the source.

Average treatment is selected with `rf_pwr=2*pi*3.0e6`, `rf_dur=1.0e-6`, and `rf_frq=-10e3`. The code calls `doublerot` with `overtone_pa` and `qnmr`, multiplies the result by `exp(1i*1.49)`, then plots its real part with `plot_1d`. The page records the inputs and output path in the example; it does not claim that the plotted trace reproduces the cited figure.
