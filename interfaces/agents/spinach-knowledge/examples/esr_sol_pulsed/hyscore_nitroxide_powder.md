# examples/esr_sol_pulsed/hyscore_nitroxide_powder.m

- Signature: `hyscore_nitroxide_powder()`

## Purpose

Time-domain Liouville-space simulation of powder-averaged HYSCORE for a ¹⁴N nitroxide at 0.350 T. The example is intended to reproduce Figure 2a of the paper cited at [doi:10.1080/00268979809483260](http://dx.doi.org/10.1080/00268979809483260). Calculation time: seconds.

## Physical model

The system contains ¹⁴N and an electron, with electron g = 2, isotropic electron–nitrogen hyperfine coupling of 5 MHz, and a nitrogen quadrupole interaction defined by `eeqq2nqi(2.4e6,0.5,1,[0 0 0])`. The simulation uses the full sphten Liouville-space basis without approximation; trajectory-level SSR and the colorbar are disabled.

## Simulation and processing

The initial state and detection operator are the electron `Lz` and `L+`, respectively. HYSCORE is calculated with `tau = 136 ns`, a 20 MHz sweep, 128 points in each dimension, and the `rep_2ang_800pts_sph` powder grid. The two-dimensional signal has its mean removed, is apodised with cosine windows in both dimensions, and is zero-filled to 256×256 before a 2D Fourier transform. The absolute spectrum is plotted with positive contours in MHz.
