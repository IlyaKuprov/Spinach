# examples/nmr_proteins/hsqc_ubiquitin_a.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_proteins/hsqc_ubiquitin_a.m)

- Signature: `hsqc_ubiquitin_a()`

## Purpose

A simulated 2D 1H-15N HSQC of human ubiquitin, with decoupling in both dimensions. The source estimates hours of calculation time and notes faster calculation with a Tesla A100 GPU. The code calculates the spectrum; it does not load an experimental HSQC spectrum or report a measured match. Source authors are Zenawi Welderufael, Luke Edwards, and Ilya Kuprov.

## Protein, spin system, and acquisition

The code calls `protein('1D3Z.pdb','1D3Z.bmrb',options)` with molecule 1, `noshift='delete'`, and `select='backbone-hsqc'`. These are protein structure and shift inputs, not a measured 2D spectrum. The field is 11.7395 T and the interaction/proximity cutoffs are 5.0/4.0 (units are not stated in the assignments).

The basis is `sphten-liouv` with IK-1, scalar-coupling connectivity, and interaction/proximity levels 4/1. The code enables `greedy`; the adjacent `gpu` text is commented out, not enabled. The sequence call is `liquid(...,@hsqc,...,'nmr')`, with spins `15N` and `1H`, and `J=90`. F1 decouples `1H` and `13C`; F2 decouples `15N` and `13C`. Sweeps are [2000, 4000], offsets [-5870, 3753], acquired points [128, 256], and zero-fill sizes [1024, 1024]. Axes are displayed in ppm.

## Processing and output

The positive and negative FIDs receive squared-cosine apodisation. After the F2 Fourier transform, the code forms the States signal as `f1_pos+conj(f1_neg)`, then performs the F1 Fourier transform. It plots the real part of the 2D spectrum.
