# examples/nmr_solids/mas_powder_ala_gridfree.m

Source: [examples/nmr_solids/mas_powder_ala_gridfree.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/mas_powder_ala_gridfree.m)

## Model

This example builds an alanine powder spin system from the PCM-DFT output in `../standard_systems/alanine.log`. The `g2spinach` call maps carbon to `13C` and nitrogen to `15N` and supplies reference shieldings `[182.1 264.5]` for zero-ppm chemical-shift references. The code sets a 14.1 T field. The individual interaction tensors are supplied by the parsed calculation output rather than listed as numbers in this script.

The source describes the target as a `13C` MAS spectrum assuming `1H` decoupling. Its experiment settings specify a rotor axis of `[1 1 1]` and a rate of 2000 Hz, and acquire on `13C` with `L+` initial and receiver operators. The script sets `parameters.decouple={}` and defines no RF pulse sequence, so the decoupling assumption in the header is not implemented as an explicit decoupling waveform here.

## Calculation and display

The basis uses `sphten-liouv` with no approximation, a `15N` longitudinal subspace, and the `+1` projection. The source sets interaction and proximity cutoffs to 5.0 and 4.0, disables trajectory-level storage, and uses maximum rank 17. It calls `gridfree(spin_system,@acquire,parameters,'nmr')`, then applies exponential apodisation with parameter 6, Fourier transforms the FID after zero-filling 256 points to 1024, and plots the real spectrum. The sweep is 50 kHz, the offset is 15 kHz, and the ppm axis is inverted. These are simulation and display settings; the file does not present an experimentally measured spectrum.
