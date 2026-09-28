# examples/nmr_liquids/coloc_test.m

- Signature: `coloc_test()`

## Purpose

A two-spin ¹H–¹³C COLOC pulse-sequence example with a long-range coupling; the source estimates seconds of calculation time.

## Physical and numerical content

The system is set to 11.7 T with shifts 4.0 and 75.0, a ¹H–¹³C scalar coupling of 5.0, and a zero self-coupling entry for spin 2. With the full sphten-Liouville basis, it simulates `liquid(...,@coloc,...,'nmr')` using delta2 = 30e-3, offsets [2250 5000], sweeps [5000 12000], and 256 × 256 points (zero-filled to 512 × 512). The FID receives cosine apodisation, then a two-dimensional Fourier transform; the plotted spectrum is its magnitude (axis units: ppm).
