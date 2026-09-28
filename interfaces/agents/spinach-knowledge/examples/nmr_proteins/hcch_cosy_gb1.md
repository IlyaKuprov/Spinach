# examples/nmr_proteins/hcch_cosy_gb1.m

- Signature: `hcch_cosy_gb1()`

## Purpose

Simulate a 3D HCCH COSY experiment on GB1 protein. The source states a calculation time of hours.

## Implementation

- Import protein data from `2N9K.pdb` and `2N9K.bmrb` with `pdb_mol=1`, `noshift='delete'`, and `select='all'`; set the magnetic field to `14.1` and the interaction and proximity cutoffs to `20.0` and `4.0`.
- Use the `sphten-liouv` basis with `IK-1` approximation, `scalar_couplings` connectivity, interaction level `4`, and proximity level `1`. Enable `greedy` and `prop_cache`.
- Set `J_ch=140`, `J_cc=35`, `delta=1.1e-3`, sweep widths `[6000 13000 6000]`, spins `{'1H','13C','1H'}`, offsets `[2500 7000 2500]`, acquisition points `[128 128 128]`, zero-fill sizes `[256 256 256]`, F3 decoupling of `13C`, and axis units of `ppm`.
- Create the spin system, remove `15N` spins, build the basis, and run `liquid(spin_system,@hcch_cosy,parameters,'nmr')`.
- Apply squared-cosine apodisation in all three dimensions to the four `pos_pos`, `pos_neg`, `neg_pos`, and `neg_neg` FIDs. Perform zero-filled, shifted FFTs along F3, F2, and F1, combining the appropriate components with complex conjugates to obtain the absorption components after F3 and F2.
- Plot `imag(spectrum)` using `plot_3d` with threshold `10`, bounds `[0.05 0.25 0.05 0.25]`, dimension `2`, and `'positive'` selection.