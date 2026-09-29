# examples/nmr_proteins/hcch_cosy_simple.m

- Signature: `hcch_cosy_simple()`
- Source: [examples/nmr_proteins/hcch_cosy_simple.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_proteins/hcch_cosy_simple.m)

## Task

A self-contained three-dimensional HCCH COSY forward simulation for a small protein-fragment model. It creates a five-spin nuclear system directly; it does not import protein files or measured spectral data, and it contains no paramagnetic centre or electron-spin interaction. The source gives a calculation time of seconds.

## Spin model and basis

The field is set to `14.1` (unit not stated). In order, the isotopes are `13C, 13C, 13C, 1H, 1H`, labelled `C, CA, CB, HA, HB`. Scalar Zeeman values are `{180 60 40 4 2}`. Scalar coupling entries are `(2,3)=35.0`, `(1,3)=0.3`, `(1,2)=55.0`, `(3,5)=130.0`, `(4,5)=7.0`, `(3,4)=6.0`, `(2,4)=140.0`, `(1,4)=4.0`, and `(5,5)=0`. The basis uses `sphten-liouv` and `approximation='none'`. Units for the field, Zeeman values, and couplings are not stated in this file.

## Sequence, acquisition, and processing

The code sets `J_ch=140`, `J_cc=35`, and `delta=1.1e-3` (units not specified). The `@hcch_cosy` sequence dimensions are `{'1H','13C','1H'}`; sweep is `[4000 7500 4000]`, offset is `[2000 7000 2000]`, acquisition points are `[64 64 64]`, and zero-fill sizes are `[256 256 256]`. It sets `decouple_f3={'13C'}` and `axis_units='ppm'`. The source does not state units for its sweep and offset literals.

The four simulated FID components `pos_pos`, `pos_neg`, `neg_pos`, and `neg_neg` are apodised with `sqcos` in all dimensions. After zero-filled shifted F3 Fourier transforms, the code forms `f3_pos=f3_pos_pos+conj(f3_neg_neg)` and `f3_neg=f3_neg_pos+conj(f3_pos_neg)`; it transforms both along F2, combines `f3f2_pos+conj(f3f2_neg)`, and transforms along F1.

## Output and scope

The plotted result is `imag(spectrum)`, with `plot_3d` threshold `10`, bounds `[0.05 0.25 0.05 0.25]`, dimension `2`, and selection `'positive'`. The source only plots the simulation; it does not save or compare an experimental spectrum. This H-C-H sequence is not an H-N-C triple-resonance experiment.
