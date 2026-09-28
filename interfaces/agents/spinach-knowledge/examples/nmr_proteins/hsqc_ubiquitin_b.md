# examples/nmr_proteins/hsqc_ubiquitin_b.m

- Signature: `hsqc_ubiquitin_b()`

## Purpose

A 1H–15N HSQC of human ubiquitin without 1H decoupling in F1 or 15N decoupling in F2, retaining nitrogen–proton multiplicity in both dimensions. The source estimates hours of calculation time and notes that a Tesla A100 GPU can make it faster. Source authors: Zenawi Welderufael, Luke Edwards, and Ilya Kuprov.

## Setup and acquisition

The example imports the backbone-HSQC selection from `1D3Z.pdb` / `1D3Z.bmrb`, at 14.1 T, with inter/proximal cutoffs 5.0/4.0. The IK-1 sphten-liouv basis uses scalar-coupling connectivity and levels 4/1; greedy optimization is enabled. It sets `J=90`, sweeps [2400, 4800], offsets [-7000, 4500], points [256, 256], zero-fill sizes [1024, 1024], and spins `15N` / `1H`. The listed decoupling is `13C` in both F1 and F2; axes are in ppm.

## Processing

The `liquid` simulation uses `@hsqc`. Positive and negative FIDs are squared-cosine apodised, transformed in F2, and combined as `f1_pos+conj(f1_neg)` for the States signal before the F1 transform. The magnitude spectrum is plotted.
