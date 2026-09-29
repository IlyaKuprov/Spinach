# examples/nmr_liquids/inept_strychnine.m

- Signature: `inept_strychnine()`
- Source: [examples/nmr_liquids/inept_strychnine.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/inept_strychnine.m)

## Purpose

An INEPT experiment on strychnine, not an INADEQUATE or NOE experiment. It models polarisation transfer between proton and carbon channels and plots the carbon response; the source estimates calculation time in minutes.

## Implementation

The spin system comes from `strychnine({'1H','13C'})`, at `5.9` T and temperature `298` K. The scalar-coupling Liouville basis uses IK-2 with proximity level 1. Sequence channels are `13C` and `1H`, with transfer coupling `J=150` Hz, sweep value `10000`, offsets `[5000 0]`, 2048 points, and zero filling to 8196 points; source units are not stated for sweep or offsets.

The code creates `13C` isotopomers and simulates each with `@inept`, exponentially apodises each FID with parameter 6, and sums their Fourier transforms. The plotted observable is the imaginary part of the summed spectrum on the `13C` channel, with the first offset used and the axis in ppm.
