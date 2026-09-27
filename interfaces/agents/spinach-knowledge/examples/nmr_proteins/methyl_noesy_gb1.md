# examples/nmr_proteins/methyl_noesy_gb1.m

- Signature: `methyl_noesy_gb1()`

## Purpose

A 1H–1H NOESY simulation of GB1 deuterated everywhere except methyl groups. Deuteria remain in the spin system because they participate in the coupling network; methyl rotation is not modelled. The source estimates hours of calculation time.

## Setup and acquisition

The example imports `2N9K.pdb` / `2N9K.bmrb`, deuterates non-methyl positions, and sets the field to 21.1356 T. It uses an inter-spin cutoff of 100 (the source comment says this retains significant dipole–dipole couplings) and a proximity cutoff of 5.0, to be increased until convergence. Redfield relaxation is selected with `tau_c=5e-9` s, `rlx_keep='kite'`, and zero equilibrium. The IK-1 sphten-liouv basis uses scalar-coupling connectivity and pairwise inter/proximal levels 2/2. After creating the system, the code removes `13C` and `15N` spins. The NOESY mixing time is 200 ms; the initial state is proton `Lz`; offset 750, sweeps [3000, 3000], points [512, 512], zero-fill sizes [2048, 2048], and axes in ppm.

## Processing

The `liquid` simulation uses `@noesy`. Cosine and sine FIDs are squared-cosine apodised, transformed in F2, and combined as `f1_cos-1i*f1_sin` for States processing before the F1 transform. The negative real 2D spectrum is plotted.
