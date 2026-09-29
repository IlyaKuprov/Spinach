# examples/nmr_liquids/dept_strychnine.m

- MATLAB implementation: [examples/nmr_liquids/dept_strychnine.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/dept_strychnine.m)

This example calculates a natural-abundance 13C DEPT135 spectrum for strychnine from a liquid-state spin model. It illustrates DEPT editing through proton-carbon scalar-coupling transfer; it does not load or fit an experimental spectrum.

## Spin system and editing

The built-in strychnine system contains 1H and 13C spins and is evaluated at 5.9 T and 298 K. The calculation uses a sphten-liouv / IK-2 basis with scalar-coupling connectivity and proximity level 1, dilutes over 13C isotopomers, then simulates each in parallel with the DEPT sequence. Its working J coupling is 150 Hz and the selection-pulse angle is 3*pi/4 radians (135 degrees). The sequence evolves for J-dependent delays, applies proton/carbon pulses and decouples the proton spins for carbon detection. The receiver state is the 13C L+ operator. At this DEPT135 angle, the sequence documentation states that CH and CH3 signals are opposite in phase to CH2 signals; quaternary carbons are not present in DEPT. The pulse-sequence reference is [DEPT](https://doi.org/10.1016/0022-2364(82)90286-4).

## Acquisition and display

The 13C sweep width is 10000 Hz with 2048 acquired points and zero filling to 8196 points. The configured offsets are [5000, 0] Hz; the plot uses the 13C offset and displays the real spectrum in ppm. Each isotopomer's FID is exponentially apodised with the source setting `{'exp',6}` before Fourier transformation and summation.

The source labels the expected calculation time as minutes. The plotted model can illustrate the DEPT135 phase edit, but the example does not supply carbon assignments, a measured spectrum, or a claimed experimental agreement.
