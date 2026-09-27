# examples/nmr_solids/case_studies/mathies_carbonate/mas_powder_mhc_fplanck.m

- Signature: `mas_powder_mhc_fplanck()`

## Purpose

Simulates the magic-angle-spinning proton NMR spectrum of all protons in the monohydrocalcite unit cell. The source cites https://doi.org/10.1038/s41467-023-44381-x and estimates hours of runtime, or minutes on a GPU.

## Model and calculation

- Reads `mhc.magres`, removes C, O, and Ca, assigns the remaining 18 sites as 1H, converts the CASTEP shielding tensors using the Huang et al. ACIE 2021 parametrisation, and uses their Cartesian coordinates.
- Uses a 9.4 T field and an `sphten-liouv` basis with approximation `IK-0`, inter-level 3, and `+1` projections. The interaction cutoff is 500 Hz. The GPU enable line is commented out.
- The pulse-acquire setup uses MAS rate 10,000 Hz, axis `[1 1 1]`, maximum rank 13, grid `rep_2ang_100pts_sph`, sweep `1/(5e-6)` Hz, 512 points, and 1024-point zero filling.
- Calls `singlerot` with `@acquire`; applies exponential apodisation (6), Fourier transforms, and plots the real spectrum.
