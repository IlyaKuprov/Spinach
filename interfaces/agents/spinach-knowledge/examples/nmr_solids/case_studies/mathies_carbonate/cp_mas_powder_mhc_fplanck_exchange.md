# examples/nmr_solids/case_studies/mathies_carbonate/cp_mas_powder_mhc_fplanck_exchange.m

- Signature: `cp_mas_powder_mhc_fplanck_exchange()`

## Purpose

Calculates cross-polarisation contact curves under magic-angle spinning for H1, H4, and C19 in monohydrocalcite, with chemical exchange between two three-spin endpoints. The source cites further details at https://doi.org/10.1038/s41467-023-44381-x and reports hours of CPU time, much faster on a GPU.

## Model and calculation

- Reads `mhc.magres`, removes O and Ca, and constructs two endpoints containing H1/H4/C19 and H4/H1/C19, respectively. The H1 and H4 sites exchange; endpoint concentrations are `[1 1]`.
- Uses the Huang et al. ACIE 2021 shift parametrisation, Cartesian coordinates, a 9.4 T field, and an `sphten-liouv` basis with no approximation. The GPU enable line is present but commented out.
- The exchange-rate series is 10, 100, 1,000, 10,000, 100,000, and 1,000,000 Hz. Each endpoint rate matrix is built from the selected rate as `[-1 1; 1 -1]`.
- Cross-polarisation settings: MAS rate 10,000 Hz, axis `[1 1 1]`, maximum rank 7, grid `rep_2ang_800pts_sph`, offsets `[2e3 1e4]` Hz, high-power field 83 kHz, and CP powers `[60 50]` kHz. Each curve uses 1,000 steps of 10 microseconds. The source calls `singlerot` with `@cp_contact_soft` and plots all six contact curves.
