# examples/nmr_liquids/deptq_strychnine.m

- MATLAB implementation: [examples/nmr_liquids/deptq_strychnine.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/deptq_strychnine.m)

This example calculates a natural-abundance 13C DEPTQ135 spectrum for strychnine. It uses the liquid-state DEPTQ variant to retain quaternary-carbon signals that are absent from the companion DEPT135 calculation; it is a simulated spin-system result, not an experimental spectrum.

## Spin system and editing

The built-in strychnine system contains 1H and 13C spins and is evaluated at 5.9 T and 298 K. A sphten-liouv / IK-2 basis with scalar-coupling connectivity and proximity level 1 is applied to each natural-abundance 13C isotopomer, which is simulated in parallel and contributes to the summed spectrum. The DEPTQ sequence uses a 150 Hz working J coupling and a beta selection-pulse angle of 3*pi/4 radians (135 degrees). Its documentation identifies this as the fixed-first-proton-pulse DEPTQ135 variant: beta controls the final proton editing pulse. The sequence evolves through J-dependent delays, applies proton and carbon pulses, decouples protons for carbon detection, and observes the 13C L+ receiver state. Unlike the DEPT sequence, its documentation explicitly notes that quaternary carbons appear. See the [DEPTQ sequence reference](https://doi.org/10.1006/jmre.1998.1595).

## Acquisition and display

The 13C sweep width is 10000 Hz, with 2048 acquired points and zero filling to 8196 points. The configured offsets are [5000, 0] Hz; the plotted axis uses the 13C offset and is in ppm. Each isotopomer's FID receives exponential apodisation with the source setting `{'exp',6}`; the real Fourier-transformed signals are summed for the one-dimensional plot.

The source labels the expected calculation time as minutes. The example supplies no carbon assignments or measured comparison, and it does not state a DEPTQ phase-sign rule for each carbon multiplicity; interpret the output as a simulated edited spectrum, not an experimental match.
