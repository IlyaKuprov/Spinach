# examples/optimal_control/features_diss_drift.m

[Stable source link](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/features_diss_drift.m) · [Background paper](http://dx.doi.org/10.1016/j.jmr.2011.07.023)

## Spin model and objective

The simulated system is a six-site protein-backbone segment, ordered 15N(Nn), 1H, 13C(CA), 13C(CB), 13C(CO), 15N(Nn+1), at 9.4 T. The source uses shifts [119.79, 8.03, 57.32, 27.71, 177.25, 115.55] ppm and scalar couplings J(1,3)=−11, J(2,3)=140, J(3,4)=35, J(3,5)=55, J(3,6)=7, and J(5,6)=−15 Hz; other entries are left uncoupled. The normalised input is 1H `Lz` (spin 2), and the target is 13C `Lz` on spin 5 (the carbonyl site). These are model parameters, not imported experimental traces.

Relaxation is included in the drift using the `t1_t2` model, zero equilibrium, and diagonal retention. The six R1 rate entries are 1.0; the R2 entries in site order are [50, 1, 1, 1, 100, 50]. The source does not annotate units for these rate arrays. A sphten-liouv `IK-0` basis with `inter_level=4` applies the source's four-spin-correlation truncation. The drift is `hamiltonian + 1i*relaxation`, after which transmitter offsets [3214, 10000, −4800] are applied for 1H, 13C, and 15N (the source gives no units for this three-value setting). Each isotope also has an offset ensemble at −100, 0, and +100 Hz.

## Pulse optimisation and test

There are three RF channels and six quadrature controls (Lx/Ly for 1H, 13C, and 15N). The 500-slice pulse uses 40 μs per slice (20 ms total), five RF levels `2π×linspace(900,1100,5)` rad/s, LBFGS GRAPE, at most 500 iterations, and the `NS` penalty with weight 0.01. Diagnostic plots are configured for correlation order, local-per-spin content, amplitudes, and spectrogram. The random initial waveform is seeded with 0.25 on the first 20 1H-y samples and the final 20 13C-y samples.

After optimisation the pulse is scaled by the mean power level and applied to the nominal spin model with `shaped_pulse_xy` (`expv-pwc`). The reported observable is the real overlap `rho_targ' * rho`. This is one simulated test propagation, not a reported ensemble-averaged score; the source prints the score at run time but records no numerical value in the file. No hardware result is claimed.
