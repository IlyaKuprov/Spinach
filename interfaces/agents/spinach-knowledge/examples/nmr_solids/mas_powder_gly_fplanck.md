# examples/nmr_solids/mas_powder_gly_fplanck.m

- MATLAB implementation: [examples/nmr_solids/mas_powder_gly_fplanck.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/mas_powder_gly_fplanck.m)

[MATLAB source](../../../../../examples/nmr_solids/mas_powder_gly_fplanck.m)

## Purpose and model

This example constructs a computed 13C MAS spectrum for glycine powder. Its header describes Fokker–Planck MAS formalism and assumes 1H decoupling; the implementation actually calls `singlerot(spin_system,@acquire,parameters,'nmr')`. The script itself does not add a 1H isotope to the imported system or set a decoupling parameter, so the header's decoupling assumption is not an explicit decoupler setting here. The source comment estimates a seconds-scale calculation time; this is not a measured runtime.

`gparse` reads `../standard_systems/glycine.log`, and `g2spinach` imports `13C` and `15N` with reference arguments `[182.1 264.5]` in that isotope order. `g2spinach` documents `references` as absolute shielding values for zero-ppm reference substances; these are calibration inputs, not reported experimental peak positions. The model field is `14.1 T`. The basis uses `sphten-liouv`, no approximation, projection `+1`, and a longitudinal `15N` subspace; interaction and proximity cutoffs are set to `5.0` and `4.0` respectively.

## MAS acquisition and spectrum

The source sets the rotor axis to `[1 1 1]`, MAS rate to `2000 Hz` (2 kHz), powder grid `leb_2ang_rank_23`, and `max_rank=23`. Acquisition uses a `5e4 Hz` sweep and 256 points, zero-fills to 1024, and sets `offset=17000`; the script does not assign `axis_units`. Both initial state and receiver are `L+` on `13C`, and the selected observed spin is `13C`.

The code calls `singlerot`, applies exponential apodisation with parameter `6`, computes `fftshift(fft(fid,parameters.zerofill))`, and plots `real(spectrum)` with `plot_1d`. The output is a simulation spectrum; this source does not supply experimentally measured output or report a numerical comparison.

Related source documentation: [g2spinach.m](../../../../../interfaces/g2spinach.m).