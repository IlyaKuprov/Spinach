# examples/nmr_solids/pdsd_simple.m

- Signature: `pdsd_simple()`
- Source: [examples/nmr_solids/pdsd_simple.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/pdsd_simple.m)

## Purpose

Calculates a two-dimensional 13C PDSD spectrum for a four-spin HCCH fragment. The source estimates minutes of calculation time and notes that a GPU is much faster.

## Spin system and rotor sampling

The model has two 13C and two 1H spins at the four coordinates listed in the source, with scalar Zeeman entries `{20.0, 5.0, 2.0, 35.0}` in isotope order and `sys.magnet=21.1356`. The source does not state units for these values. It uses the full spherical-tensor Liouville basis (`sphten-liouv`, `approximation='none'`). The simulation calls `singlerot` with the PDSD callback, rate 10000, rotor axis `[sqrt(2/3), 0, sqrt(1/3)]`, spherical grid `rep_2ang_100pts_sph`, and maximum rank 11. No gradient is configured.

## Acquisition and processing

The detected spins are 1H and 13C; the mixing-time parameter is `10e-3`. The two dimensions use 256 points each and zero-fill to 1024 each; sweep is 10000, and offsets are `[3150, 6750]`. The source does not annotate units for these acquisition values. `singlerot` returns cosine and sine signals; each receives squared-cosine apodisation in both dimensions. The code Fourier-transforms the indirect-dimension signals, forms the States combination `f1_cos - 1i*f1_sin`, transforms the direct dimension, and takes the real part for the plotted spectrum. The plotting parameters select 13C and positive contours.

Zero track elimination is explicitly enabled with `zte` in `sys.enable`.
