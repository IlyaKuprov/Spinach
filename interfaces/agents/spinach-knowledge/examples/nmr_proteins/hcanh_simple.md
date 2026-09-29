# examples/nmr_proteins/hcanh_simple.m

- Signature: `hcanh_simple()`
- Source: [examples/nmr_proteins/hcanh_simple.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_proteins/hcanh_simple.m)

## Task

A self-contained forward simulation of a three-dimensional H(CA)NH protein-NMR-style sequence on a five-spin toy system. This source does not import a protein structure, measured spectrum, or paramagnetic centre; it constructs a nuclear spin model directly.

## Spin model and basis

The field is set to `sys.magnet=14.1` (the source gives no unit). The spins, in order, are `15N, 13C, 1H, 13C, 1H`, labelled `N, CA, H, C, HA`. The scalar Zeeman values are `{110 60 8 180 4}`. Scalar coupling entries by one-based spin pair are `(1,3)=92`, `(1,2)=11`, `(1,4)=15`, `(1,5)=1`, `(2,4)=55`, `(2,3)=2`, `(2,5)=140`, `(3,4)=4`, `(3,5)=8`, and `(4,5)=4`. The basis uses `sphten-liouv` with `approximation='none'`. Units for the field, Zeeman values, and couplings are not stated in this file.

## Sequence, acquisition, and processing

The sequence is `@hcanh` with dimension spins `{'1H','15N','1H'}`. The source sets sweeps to `[5000 3000 5000]`, offsets to `[3600 -6600 3600]`, acquisition points to `[64 64 64]`, and zero-fill sizes to `[256 256 256]`; it labels the axes `ppm`. Units for the sweep and offset literals are not specified in the source.

The call `liquid(spin_system,@hcanh,parameters,'nmr')` generates simulated FIDs. Each of `pos_pos`, `pos_neg`, `neg_pos`, and `neg_neg` is apodised with `sqcos` in all three dimensions. The code zero-fills and Fourier-transforms along F3, forms `f3_pos=f3_pos_pos+conj(f3_neg_neg)` and `f3_neg=f3_neg_pos+conj(f3_pos_neg)`, transforms along F2 and combines `f3f2_pos+conj(f3f2_neg)`, then transforms along F1.

## Output and scope

The script displays `real(spectrum)` with `plot_3d`, threshold `10`, bounds `[0.05 0.5 0.05 0.5]`, dimension `2`, and selection `'positive'`; it contains no spectrum-file export. The source describes calculation time as seconds. No measured-data agreement is claimed.
