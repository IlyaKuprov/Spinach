# examples/nmr_liquids/cosy45_rotenone.m

- Signature: `cosy45_rotenone()`

## Purpose

Simulates a magnitude-mode 45° COSY spectrum for a 22-proton rotenone spin system. The source estimates minutes of calculation time.

## Physical and numerical content

The source specifies the proton shifts and scalar couplings, and uses an IK-2 scalar-coupling basis with proximity level 1, cutoff 4.0, and three S3 symmetry groups on spins [14–16], [17–19], and [20–22]. At 5.9 T it sets the COSY angle to pi/4, offset 1200, sweep 2000, and 512 × 512 acquired points, zero-filled to 2048 × 2048. After cosine apodisation, it computes a two-dimensional Fourier transform and plots abs(spectrum) (axis units: ppm).

Source citation: [DOI 10.1002/jhet.5570250160](http://dx.doi.org/10.1002/jhet.5570250160).
