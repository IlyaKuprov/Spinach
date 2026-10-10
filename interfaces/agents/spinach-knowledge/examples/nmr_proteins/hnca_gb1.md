# examples/nmr_proteins/hnca_gb1.m

- Signature: `hnca_gb1()`
- Source: [examples/nmr_proteins/hnca_gb1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_proteins/hnca_gb1.m)

## Task and input

A forward simulation of a three-dimensional HNCA spectrum for GB1. The source assumes only the backbone is `13C,15N`-labelled. It calls `protein('2N9K.pdb','2N9K.bmrb',options)` with `pdb_mol=1`, `noshift='delete'`, and `select='backbone-minimal'`. These inputs construct the protein spin system; the script does not load a measured spectrum or compare the simulated result with experimental data. No paramagnetic centre or electron-spin interaction is specified.

## Spin model and basis

The field literal is `14.1` (unit not stated). The interaction and proximity cutoffs are `2.0` and `4.0`. The basis uses `sphten-liouv`, approximation `IK-1`, connectivity `scalar_couplings`, interaction level `4`, and proximity level `1`. The active algorithmic options enable `zte` and `greedy` and disable `krylov`; the source's `% 'gpu'` text is a comment, not an enabled option.

## Sequence, acquisition, and processing

The `@hnca` sequence dimensions are `{'15N','13C','1H'}`. The source sets sweep `[2800 5000 3000]`, offset `[-7200 8600 5100]`, acquisition points `[128 128 128]`, zero-fill sizes `[256 256 256]`, and `axis_units='ppm'`. The source does not state units for the sweep or offset literals.

The `liquid(spin_system,@hnca,parameters,'nmr')` call generates simulated FIDs. Each of `pos_pos`, `pos_neg`, `neg_pos`, and `neg_neg` is apodised with `sqcos` in all dimensions. The zero-filled, shifted F3 transforms are combined as `f3_pos=f3_pos_pos+conj(f3_neg_neg)` and `f3_neg=f3_neg_pos+conj(f3_pos_neg)`; after the F2 transforms the source combines `f3f2_pos+conj(f3f2_neg)`, then transforms along F1.

## Output and scope

The plot is `-real(spectrum)`, using threshold `10`, bounds `[0.1 0.5 0.1 0.5]`, dimension `2`, and selection `'positive'`. The source says calculation takes minutes and is faster with a Tesla A100 GPU; it does not activate the commented `gpu` option. No spectrum-file export is coded.
