# examples/dnp_liq/ccdnp/rates_main_text.m

- MATLAB implementation: [examples/dnp_liq/ccdnp/rates_main_text.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_liq/ccdnp/rates_main_text.m)

## What it computes

Builds the Redfield relaxation superoperator for a liquid-state cross-correlated DNP model with one proton and two exchange-coupled electrons, then prints selected self- and cross-relaxation projections. Use it to inspect the listed operator-specific rates for this main-text parameter set; it is not a frequency-scan or polarisation time-course calculation.

## Running and inputs

Call the no-argument MATLAB function `rates_main_text()` with Spinach on the MATLAB path. It constructs the spin system, basis and relaxation superoperator internally; it accepts no parameters and reads no external input files. The source comment estimates a calculation time of seconds. It uses a full `sphten-liouv` basis (`bas.approximation='none'`) and Redfield relaxation.

## Model settings

- Isotopes: `{'1H','E','E'}`. The file sets `sys.magnet=14.1` T. Proton Zeeman chemical shifts are `[0 10 20]` ppm with Euler triple `[0 0 0]`.
- Electron Zeeman principal values are `[2.003400 2.003800 2.003800]` and `[2.005700 2.003000 2.003000]`, with Euler triples `[-0.872 -0.013 0.868]` and `[-1.145 0.061 1.143]`. Scalar electron–electron coupling is `3e6` Hz.
- Coordinates are `[0 0 0]`, `[5.090 0.010 0.958]`, and `[-5.090 0.061 1.032]`. Coordinates are in ångström (Å).
- Relaxation settings: `inter.relaxation={'redfield'}`, `inter.equilibrium='zero'`, `inter.rlx_keep='labframe'`, `inter.temperature=298`, `inter.tau_c={100e-12}`, and `sys.tols.rlx_integration=1e-10`. Temperature is in kelvin and `tau_c` in seconds; `rlx_integration` is a numerical tolerance.

## Console output

After computing `R=relaxation(spin_system)`, the script prints matrix-element projections for normalised longitudinal proton/electron operators (their self terms and electron-to-proton terms), four normalised electron transverse combinations of the form `E1p ± 2*E1pE2z` and `E2p ± 2*E1zE2p`, and proton-containing longitudinal cross terms (including the printed terms labelled for NzE1z, NzE2z, and NzE1zE2z to Nz). It also prints mixed transverse-coherence projections labelled E1p to NzE1p, E1p to NzE1pE2z, and E1pE2z to NzE1pE2z. The source-specific clarification is that these are selected contractions of the relaxation superoperator with constructed operator states—not a dump of the complete relaxation matrix or its eigenmodes. It writes no data file and creates no figure; the values appear only in MATLAB's command output.

## Reference and caveat

The source cites https://doi.org/10.1016/j.jmr.2021.106940. numerical output depends on running the MATLAB/Spinach calculation.
