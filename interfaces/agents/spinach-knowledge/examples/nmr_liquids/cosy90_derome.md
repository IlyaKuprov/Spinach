# examples/nmr_liquids/cosy90_derome.m

- Signature: `cosy90_derome()`

## Purpose

Reproduces Figure 8.26 from Andrew Derome's *Modern NMR Techniques for Chemistry Research* with a three-proton 90° COSY simulation. The source estimates seconds of calculation time.

## Physical and numerical content

The three proton shifts are 3.70, 3.92, and 4.50; the specified scalar couplings are J12 = 10, J23 = 12, and J13 = 4. At 16.1 T, the source uses the full sphten-Liouville basis and sets the COSY angle to pi/2, offset 2800, sweep 700, and 1024 × 1024 points with 2048 × 2048 zero filling. The simulated FID is square-cosine apodised, Fourier transformed in both dimensions, and the real spectrum is plotted (axis units: ppm).
