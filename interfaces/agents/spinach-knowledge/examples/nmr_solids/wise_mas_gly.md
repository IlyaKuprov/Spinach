# examples/nmr_solids/wise_mas_gly.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/wise_mas_gly.m)

- Signature: wise_mas_gly()

## Purpose

Simulates a two-dimensional WISE spectrum of alpha-glycine powder under magic-angle spinning (MAS). The source estimates a runtime of hours and notes faster execution on a GPU.

## Spin system and rotor sampling

The seven-spin system is generated from the PCM-DFT glycine structure: two 13C and five 1H sites. The file labels the setting as a 400 MHz spectrometer and sets sys.magnet to 9.4. Isotropic shifts are set to 176.4 and 43.6 for the carbonyl and Cα carbons, 2.6 and 3.8 for the two Hα sites, and 8.0 for each of three HN sites. It uses an sphten-liouv basis with IK-0 approximation and inter-level 3, discarding interactions below 200 Hz. The rotor-axis vector is [1, 1, 1], the rate parameter is 5000, and the orientation grid is rep_2ang_200pts_sph; maximum rank is 9.

## WISE acquisition and processing

The wise sequence uses 1H and 13C channels, offsets [2000, 10000] Hz, high-power irradiation of 83000 Hz, cross-polarisation powers [60000, 50000] Hz, and a 0.0001-second contact duration. The two sweeps are [1/(6e-6), 1/(33e-6)] Hz, with [128, 512] acquired points and [512, 2048] zero-filled points. Detection is on the 13C L+ state.

singlerot returns cosine and sine FIDs. Each is apodised with a squared-cosine window; the first-dimension transforms are combined as real cosine plus i times real sine, then Fourier transformed in the second dimension. The plotted observable is the real part of the resulting two-dimensional spectrum.

Zero track elimination is explicitly enabled with `zte` in `sys.enable`.
