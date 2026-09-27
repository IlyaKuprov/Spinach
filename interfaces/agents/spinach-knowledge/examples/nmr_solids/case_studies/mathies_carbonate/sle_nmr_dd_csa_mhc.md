# examples/nmr_solids/case_studies/mathies_carbonate/sle_nmr_dd_csa_mhc.m

- Signature: `sle_nmr_dd_csa_mhc()`

## Purpose

Calculates proton NMR spectra for water protons in monohydrocalcite under MAS while varying slow isotropic rotational-diffusion correlation time. The source cites https://doi.org/10.1038/s41467-023-44381-x and reports minutes of runtime, or seconds with a GPU.

## Model and calculation

- Reads `mhc.magres`, removes C, O, and Ca, and retains the proton sites at positions 1 and 4. The CASTEP shielding tensors are converted to shifts using the Huang et al. ACIE 2021 parametrisation; the model uses their coordinates, a 9.4 T field, and an `sphten-liouv` basis with no approximation.
- Uses `gridfree` acquisition with MAS rate 10,000 Hz, axis `[1 1 1]`, sweep 120,000 Hz, 1,024 points, 4,096-point zero filling, zero offset, and no decoupling.
- Runs correlation times `1e-6*[0.10 1.00 10.0 100.0 1000.0]` s with corresponding Wigner ranks `[2 3 5 7 13]`. For each, it acquires the FID, applies exponential apodisation (6), Fourier transforms, and plots the real spectrum. GPU enablement is commented out in the source.
