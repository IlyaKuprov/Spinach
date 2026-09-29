# examples/nmr_proteins/hcch_cosy_gb1.m

- Signature: `hcch_cosy_gb1()`
- Source: [examples/nmr_proteins/hcch_cosy_gb1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_proteins/hcch_cosy_gb1.m)

## Task and input

A three-dimensional HCCH COSY forward simulation using GB1 protein input. The code calls `protein('2N9K.pdb','2N9K.bmrb',options)` with `pdb_mol=1`, `noshift='delete'`, and `select='all'`; it then removes all imported `15N` spins with `kill_spin`. These are structure/chemical-shift inputs to construct the spin system, not imported measured FIDs or a measured spectrum. The source does not encode a paramagnetic centre or electron-spin interaction.

## Spin model and basis

The field literal is `14.1` (unit not stated). Interaction and proximity cutoffs are `20.0` and `4.0`. The basis is `sphten-liouv`, approximation `IK-1`, connectivity `scalar_couplings`, interaction level `4`, and proximity level `1`. Algorithmic options enabled are `greedy` and `prop_cache`.

## Sequence, acquisition, and processing

The code sets `J_ch=140`, `J_cc=35`, and `delta=1.1e-3`; it does not state units for these literals. The three sequence dimensions use `{'1H','13C','1H'}`, with sweep `[6000 13000 6000]`, offset `[2500 7000 2500]`, acquisition points `[128 128 128]`, and zero-fill sizes `[256 256 256]`. It sets `decouple_f3={'13C'}` and `axis_units='ppm'`; units for sweep and offset are not stated in the source.

The call `liquid(spin_system,@hcch_cosy,parameters,'nmr')` produces simulated FIDs. The four components `pos_pos`, `pos_neg`, `neg_pos`, and `neg_neg` are each apodised with `sqcos` in all dimensions. The zero-filled, shifted F3 transforms are combined as `f3_pos=f3_pos_pos+conj(f3_neg_neg)` and `f3_neg=f3_neg_pos+conj(f3_pos_neg)`; the F2 transforms are combined as `f3f2=f3f2_pos+conj(f3f2_neg)`, followed by the F1 transform.

## Output and scope

The script displays `imag(spectrum)` with `plot_3d`, threshold `10`, bounds `[0.05 0.25 0.05 0.25]`, dimension `2`, and selection `'positive'`; no spectrum-file export is coded. The source estimates calculation time in hours. This is an H-C-H three-dimensional experiment, not an H-N-C triple-resonance sequence; no experimental spectrum comparison is coded.
