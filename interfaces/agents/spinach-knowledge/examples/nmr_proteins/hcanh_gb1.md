# examples/nmr_proteins/hcanh_gb1.m

Source: [examples/nmr_proteins/hcanh_gb1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_proteins/hcanh_gb1.m)

## Purpose

Simulates a three-dimensional H(CA)NH protein NMR spectrum for GB1, with the source comment's assumption that only the backbone is 13C,15N-labelled. It is a Spinach simulation, not processing of measured data. The sequence source identifies F1 as 1H, F2 as 15N, and F3 as 1H. No paramagnetic centre or magnetic tensor is specified in these example or sequence sources; tensor units and experimental agreement therefore do not arise here.

## Protein system and basis

The example imports `2N9K.pdb` and `2N9K.bmrb` using `pdb_mol=1`, `noshift='delete'`, and `select='backbone-minimal'`, sets the field to `14.1` T, and sets `inter_cutoff=2.0` and `prox_cutoff=4.0` (the source does not state units for these tolerance values). It builds a `sphten-liouv` basis with `IK-1` approximation, `scalar_couplings` connectivity, `inter_level=4`, and `prox_level=1`. `sys.enable={'greedy'}` is active; `gpu` appears only as a commented alternative.

## Sequence, acquisition and processing

The example calls `liquid(spin_system,@hcanh,parameters,'nmr')` with spins `1H`, `15N`, `1H`, sweeps `[6000 3000 6000]` Hz, offsets `[4200 -7200 4200]` Hz, 128 acquired points per dimension, and zero filling `[256 256 256]`; display axes are ppm. In `experiments/nmr_protein/hcanh.m`, the three dimensions are F1=`1H`, F2=`15N`, F3=`1H`. Its hard-coded transfer couplings are `J_ch=140` Hz and `J_nh=92` Hz; `delta2=12.5e-3` s and `delta4=23.0e-3` s, while other listed transfer delays are computed as `1/(4J)` in seconds. The code selects positive and negative 1H coherence for F1 and positive and negative 15N coherence for F2, returning four States-channel FIDs: `pos_pos`, `pos_neg`, `neg_pos`, and `neg_neg`.

All four FIDs receive squared-cosine apodisation. The example Fourier-transforms F3 for each channel; for F2 it transforms the positive and negative channel pairs separately and combines them with complex conjugation as `f3f2_pos + conj(f3f2_neg)`; it then Fourier-transforms F1. The real part is displayed with `plot_3d`; no output data file is written.

The sequence source cites the bidirectional-propagation method at http://dx.doi.org/10.1016/j.jmr.2014.04.002 and documents the H(CA)NH sequence as Figure 7.37 in *Protein NMR Spectroscopy*, 2nd edition. See also [experiments/nmr_protein/hcanh.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/nmr_protein/hcanh.m).
