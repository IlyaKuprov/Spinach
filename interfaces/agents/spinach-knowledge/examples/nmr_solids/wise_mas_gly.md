# examples/nmr_solids/wise_mas_gly.m

- Signature: `wise_mas_gly()`

## Purpose

Simulates a two-dimensional WISE (wide-line separation) spectrum of alpha-glycine powder under magic-angle spinning. The source notes an hours-long runtime, substantially shorter on a GPU.

## Spin system and experiment

The PCM-DFT spin system is loaded from `../../examples/standard_systems/glycine.log`; the field is 9.4 T. Isotropic shifts are set to 176.4 and 43.6 ppm for CO and Cα, 2.6 and 3.8 ppm for the two Hα sites, and 8.0 ppm for each of the three HN sites. The calculation uses the `sphten-liouv` basis with IK-0 approximation and inter-level 3. MAS rate is 5 kHz, the axis is [1 1 1], maximum rank is 9, and the powder grid is `rep_2ang_200pts_sph`.

The WISE/CP settings are offsets [2000, 10000] Hz, high-power irradiation 83 kHz, CP powers [60, 50] kHz, and CP duration 100 μs. The two dimensions use sweeps [1/(6 μs), 1/(33 μs)] Hz, [128, 512] acquired points, and [512, 2048] zero-filled points; the listed spin channels are 1H and 13C.

## Processing

The sequence is simulated by `singlerot` with the WISE sequence function. Cosine and sine components are Fourier transformed along the first dimension and combined as `real(F1_cos) + i real(F1_sin)`; a second Fourier transform produces the plotted 2D spectrum.
