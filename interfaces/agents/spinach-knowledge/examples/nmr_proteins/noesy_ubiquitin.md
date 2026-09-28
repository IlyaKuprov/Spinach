# examples/nmr_proteins/noesy_ubiquitin.m

- Signature: `noesy_ubiquitin()`

## Purpose

A 1H–1H NOESY spectrum of ubiquitin with a 65 ms mixing time; the protein is assumed not to be 13C- or 15N-labelled. The source estimates hours of calculation time.

## Setup and acquisition

The example imports `1D3Z.pdb` / `1D3Z.bmrb` with all atoms selected, at 21.1356 T. It uses inter/proximal cutoffs 2.0/4.0, Redfield relaxation with `tau_c=5e-9` s, `rlx_keep='kite'`, zero equilibrium, and an IK-1 sphten-liouv basis with scalar-coupling connectivity and levels 4/3. Propagation caching and greedy optimization are enabled; the code removes `13C` and `15N` spins. The acquisition uses `tmix=0.065` s, proton `Lz` initial state, offset 4250, sweeps [11750, 11750], points [512, 512], zero-fill sizes [2048, 2048], and ppm axes.

## Processing

The `liquid` simulation uses `@noesy`. Cosine and sine FIDs are squared-cosine apodised, transformed in F2, and combined as `f1_cos-1i*f1_sin` for States processing before the F1 transform. The negative real 2D spectrum is plotted.
