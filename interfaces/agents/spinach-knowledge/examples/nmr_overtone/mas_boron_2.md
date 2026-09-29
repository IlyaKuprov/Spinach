# examples/nmr_overtone/mas_boron_2.m

- MATLAB implementation: [examples/nmr_overtone/mas_boron_2.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_overtone/mas_boron_2.m)

This example simulates a 10B MAS overtone spectrum with the sample spinning in the JEOL direction. The source credits the parameters to Nghia Duong and Yusuke Nishiyama and describes the RF power and pulse width as realistic; it estimates hours of calculation. The file supplies no DOI or measured comparison, so those comments are not evidence of experimental agreement.

The model uses isotope 10B, magnet setting 16.4, and the coupling input `eeqq2nqi(0.7e6,0.0,3,[0 0 0])`. Units are not annotated for the magnet setting or coupling arguments. Relaxation is `damp` with diagonal retention, zero equilibrium, and `damp_rate=100`; the basis is `sphten-liouv` with approximation `none`.

The MAS settings are rank 12, axis `[sqrt(2/3) 0 sqrt(1/3)]`, rate 70000, and grid `rep_2ang_200pts_sph`. The spectrum uses sweep `[-141e3 -139e3]`, 256 points, 256-point zero-fill, and `axis_units='kHz'`. The initial state is 10B `Lz`. The source defines the receiver as `cos(theta)*Lz state + sin(theta)*Lx state`, and defines a corresponding Lx operator using the same weights, with `theta=atan(sqrt(2))`.

Average treatment uses `rf_pwr=2*pi*50e3/sin(theta)`, `rf_dur=2e-3`, and `rf_frq=-140e3`. The code calls `singlerot` with `overtone_pa` and `qnmr`, multiplies the spectrum by `exp(1i*1.45)`, and plots its real part using `plot_1d`. The source values are recorded as written; only the spectral-axis setting is explicitly labelled kHz.
