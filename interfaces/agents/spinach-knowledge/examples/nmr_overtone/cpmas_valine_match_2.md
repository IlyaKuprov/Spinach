# examples/nmr_overtone/cpmas_valine_match_2.m

- Signature: `cpmas_valine_match_2()`

## Purpose

Maps the Hartmann–Hahn condition for proton-to-`14N` overtone cross-polarisation in N-acetylvaline under MAS, using the Fokker–Planck formalism. The calculation scans both proton RF power and sample spinning rate over a rough powder grid; the source estimates hours of calculation time.

The valine quadrupolar tensor data are attributed to [the cited paper](https://doi.org/10.1039/c4cp03994g).

## Model and calculation

The system uses `14N` and `1H` at 14.10220742 T, with nitrogen quadrupole parameters 3.21 MHz and asymmetry 0.27 (spin 1), shifts [57.5, 81.0, 227.0] for nitrogen and zero for protons. Diagonal damping relaxation has rate 10000; the basis is `sphten-liouv` without approximation.

The average-treatment calculation uses rank 5, grid `rep_2ang_200pts_oct`, 256 points and 256-point zero-fill, and sweeps ±15 kHz about the RF frequency. It scans 50 proton-power values from 10 to 200 kHz and 50 spinning rates from 20 to 90 kHz. For each pair, the code sets the proton RF field using 55 kHz and the selected power, sets the spinning rate to the negative of the sampled rate, and sets the RF frequency to 46.30 kHz minus twice that signed rate. It runs `singlerot` with `overtone_cp` in a `parfor` loop and records the summed real spectrum intensity.
