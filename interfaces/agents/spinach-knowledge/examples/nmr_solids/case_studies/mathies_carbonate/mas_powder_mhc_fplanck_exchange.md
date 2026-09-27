# examples/nmr_solids/case_studies/mathies_carbonate/mas_powder_mhc_fplanck_exchange.m

- Signature: `mas_powder_mhc_fplanck_exchange()`

## Purpose

Simulates MAS proton NMR for water protons in monohydrocalcite with position exchange between two reaction endpoints. The source cites https://doi.org/10.1038/s41467-023-44381-x and reports seconds of runtime.

## Model and calculation

- Reads `mhc.magres`, removes C, O, and Ca, and forms two endpoints from the proton sites at positions 1 and 4, swapped between endpoints. Their concentrations are `[1 1]`; the two-state exchange-rate matrix uses 2,000 Hz.
- Uses a 9.4 T field, the Huang et al. ACIE 2021 shielding-to-shift parametrisation, an `sphten-liouv` basis with no approximation, and the selected Cartesian coordinates.
- Acquisition settings are MAS rate 10,000 Hz, axis `[1 1 1]`, maximum rank 13, grid `rep_2ang_100pts_sph`, sweep `1/(5e-6)` Hz, 512 points, and 1024-point zero filling. The source uses `singlerot` with `@acquire`, exponential apodisation (6), and a Fourier transform.
