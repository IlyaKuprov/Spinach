# examples/nmr_proteins/hncaco_gb1.m

- Signature: `hncaco_gb1()`

## Purpose

Simulates an HNCACO spectrum of GB1, assuming only the backbone is 13C,15N-labelled. The source estimates a calculation time of minutes, faster with a Tesla A100 GPU.

## Model and simulation

- Imports protein data from `2N9K.pdb` and `2N9K.bmrb` using `protein`, with `pdb_mol=1`, `noshift='delete'`, and `select='backbone-minimal'`.
- Sets the magnetic field to `14.1` and the interaction and proximity cutoffs to `2.0` and `4.0`.
- Uses the `sphten-liouv` formalism with `IK-1` approximation, `scalar_couplings` connectivity, interaction level `4`, and proximity level `1`. The `greedy` option is enabled; **Krylov is disabled**. The GPU option appears only in a comment, not in the enabled options.
- Creates the spin system and basis, then runs `liquid(spin_system,@hncaco,parameters,'nmr')`. Sequence parameters are `J_nh=92`, `T=25e-3`, `delta2=3e-3`, spins `{'15N','13C','1H'}`, sweep widths `[3000 2500 3000]`, offsets `[-7200 26500 5100]`, acquisition points `[64 64 64]`, zero-filled points `[256 256 256]`, and axis units `ppm`.

## Processing and output

Each of the four pathway signals (`pos_pos`, `pos_neg`, `neg_pos`, `neg_neg`) receives squared-cosine apodisation in all three dimensions. The code applies zero-filled, shifted FFTs successively along F3, F2, and F1. Between transforms, it combines conjugate pathway signals to obtain the absorption components: `f3_pos=f3_pos_pos-conj(f3_neg_neg)`, `f3_neg=f3_neg_pos-conj(f3_pos_neg)`, and `f3f2=f3f2_pos-conj(f3f2_neg)`. Finally, it plots `imag(spectrum)` with `plot_3d(...,10,[0.2 0.9 0.2 0.9],2,'positive')`.