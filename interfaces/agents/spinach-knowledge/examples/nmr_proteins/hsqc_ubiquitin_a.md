# examples/nmr_proteins/hsqc_ubiquitin_a.m

- Signature: `hsqc_ubiquitin_a()`

## Purpose

A 1H–15N HSQC of human ubiquitin with decoupling in both dimensions. The source estimates hours of calculation time and notes that a Tesla A100 GPU can make it faster. Source authors: Zenawi Welderufael, Luke Edwards, and Ilya Kuprov.

## Setup and acquisition

The example imports the backbone-HSQC selection from `1D3Z.pdb` / `1D3Z.bmrb`, at 11.7395 T, with inter/proximal cutoffs 5.0/4.0. The IK-1 sphten-liouv basis uses scalar-coupling connectivity and levels 4/1; greedy optimization is enabled. It sets `J=90`, sweeps [2000, 4000], offsets [-5870, 3753], points [128, 256], zero-fill sizes [1024, 1024], and spins `15N` / `1H`. F1 decouples `1H` and `13C`; F2 decouples `15N` and `13C`; axes are in ppm.

## Processing

The `liquid` simulation uses `@hsqc`. The positive and negative FIDs receive squared-cosine apodisation; after the F2 transform the States signal is formed as `f1_pos+conj(f1_neg)`, then Fourier-transformed in F1. The real 2D spectrum is plotted.
