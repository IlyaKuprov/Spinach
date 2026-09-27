# examples/nmr_solids/fitting/vanadium_csa_nqi/vfit_4sim_ik.m

- Signature: `vfit_4sim_ik()`

## Purpose

Simultaneous fitting of multiple 51V MAS NMR spectra with respect to the chemical shielding anisotropy and quadrupole coupling tensor parameters. Calculation time: hours, much faster with a GPU.

## Physical / mathematical content

The model is a single 51V spin with a chemical-shielding tensor and a quadrupolar coupling tensor. One common set of tensor parameters is used to fit four experimental spectra acquired at different MAS rates; the source lists the rates as 41, 38.5, 36, and 34 kHz.

## Numerical / algorithmic content

The script loads and Savitzky–Golay filters `v12_29_dec15.spc`, `v12_31_dec15.spc`, `v12_33_dec15.spc`, and `v12_35_dec15.spc`, extracts and normalises the selected spectral ranges, then minimises a summed squared residual with `fminsearch`. Each trial simulates the four spectra using the same 51V spin system, a rank-30 truncation, and `rep_2ang_200pts_oct`; each signal is Gaussian-apodised before Fourier transformation.

## Implementation structure

Preprocesses four experimental spectra, maps the fitted parameters to the chemical-shift and quadrupolar tensors, performs four MAS simulations at their respective rates, and compares the resulting spectra with the measurements in a four-panel plot.
