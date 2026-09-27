# examples/nmr_overtone/cpmas_valine_match_1.m

- Signature: `cpmas_valine_match_1()`

## Purpose

Simulates proton-to-`14N` overtone cross-polarisation in N-acetylvaline under magic-angle spinning (MAS), using the Fokker–Planck formalism. It calculates a Hartmann–Hahn condition profile over a rough powder grid as a function of proton RF power; the source estimates hours of calculation time.

The valine quadrupolar tensor data are attributed to [the cited paper](https://doi.org/10.1039/c4cp03994g).

## Model and calculation

The system uses `14N` and `1H` at 14.10220742 T. The nitrogen quadrupole parameters are 3.21 MHz, asymmetry 0.27, and spin 1; the listed nitrogen shifts are [57.5, 81.0, 227.0], with the proton shifts set to zero. The model uses diagonal damping relaxation (rate 2000) and the `sphten-liouv` basis without approximation.

At the magic angle, the simulation uses rank 9, a 19.840 kHz spinning rate, grid `rep_2ang_800pts_sph`, and a 70–105 kHz sweep with 256 points and 256-point zero-fill. The `14N` overtone experiment uses average treatment, 86.30 kHz RF frequency, and a 70 μs RF duration. Fifteen proton-power settings span 26–40 kHz; the proton RF field is set from each value together with 55 kHz. The spectrum is calculated with `singlerot` and `overtone_cp`.
