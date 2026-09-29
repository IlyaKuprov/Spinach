# examples/nmr_solids/mas_powder_gly_gridfree.m

- MATLAB implementation: [examples/nmr_solids/mas_powder_gly_gridfree.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/mas_powder_gly_gridfree.m)

[MATLAB source](../../../../../examples/nmr_solids/mas_powder_gly_gridfree.m)

## Purpose and model

This example constructs a computed 13C MAS spectrum for glycine powder with the grid-free Fokker–Planck route: the source calls `gridfree(spin_system,@acquire,parameters,'nmr')`. Its header says the magnetic parameters are estimated from a DFT calculation and estimates minutes for calculation time; that duration is a source estimate, not a measured runtime.

`gparse` reads `../standard_systems/glycine.log`; `g2spinach` imports `13C` and `15N` with reference arguments `[182.1 264.5]` in isotope order. The `g2spinach` argument `references` denotes absolute shielding values for zero-ppm reference substances, not measured spectrum peaks. The source sets field `14.1 T`. The basis is `sphten-liouv`, approximation `none`, longitudinal `15N`, and projection `+1`; interaction and proximity cutoffs are `5.0` and `4.0`.

## MAS acquisition and spectrum

The MAS rotor axis is `[1 1 1]` and rate is `2000 Hz` (2 kHz); powder settings are `leb_2ang_rank_23` and `max_rank=23`. The acquisition sweep is `5e4 Hz` with 256 points, zero-fill 1024, and `offset=17000`. It selects `13C`, assigns `axis_units='ppm'`, sets `invert_axis=1`, and explicitly leaves `decouple={}`. Both initial state and receiver are `L+` on `13C`.

After `gridfree` returns the FID, the script applies exponential apodisation parameter `6`, Fourier transforms with `fftshift(fft(fid,parameters.zerofill))`, and plots the real spectrum via `plot_1d`. These are model and processing settings for a computed spectrum; this file does not provide experimental measured output or a numerical comparison.

Related source documentation: [g2spinach.m](../../../../../interfaces/g2spinach.m).