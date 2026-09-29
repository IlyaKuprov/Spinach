# examples/relaxation_theory/trosy_methyl.m

- Signature: `trosy_methyl()`

## Purpose and model

Simulates methyl TROSY for a rapidly rotating `13CH3` group in a slowly tumbling protein using the source's Fokker-Planck formalism. The source estimates a calculation time of minutes and sets the field parameter to `14.1`; no unit is annotated. It represents the three methyl orientations as three four-spin rotamers, for a 12-spin system containing one carbon and three protons in each rotamer. The proton positions are cyclically permuted between rotamers, their populations are equal, and the source supplies a three-state jump-rate matrix. It sets `tau_m = 1e-11` and `k_jump = 1/(2*tau_m)`; no unit is annotated for these values.

The shielding tensors are source-provided DFT values converted to chemical-shift tensors. The three methyl-proton shifts are described in the source as guesses and adjusted by `0.8`, `1.0`, and `1.2`. Intramethyl scalar-coupling entries are set to `125` for each carbon-proton pair and `-12` for each proton-proton pair; the source does not annotate units for these entries. The system uses the full `sphten-liouv` basis, disables `zte`, and sets `tau_c = 50e-9` (no unit is annotated for this parameter) with maximum rank 3.

## Simulated spectra

Frequency-domain `gridfree` detection with `slowpass` produces separate carbon and proton spectra. Each channel uses its own `L+` initial state and matching receiver, with no decoupling. The carbon spectrum uses a sweep from -300 to 300 Hz and 1024 points; the proton spectrum uses 200 to 1000 Hz and 2048 points. The plotted outputs are the real parts of the calculated spectra, not measured spectra.

Source: [examples/relaxation_theory/trosy_methyl.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/trosy_methyl.m).
