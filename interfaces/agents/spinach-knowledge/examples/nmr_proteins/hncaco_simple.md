# examples/nmr_proteins/hncaco_simple.m

- Signature: `hncaco_simple()`

## Purpose

A minimal example of HNCACO pulse sequence simulation. Calculation time: seconds.

## Spin system and simulation

- Sets the magnet field to 14.1 and defines four spins: `15N` (N), `13C` (CA), `1H` (H), and `13C` (C).
- Sets scalar Zeeman values to `[110 60 7 180]`. Scalar couplings are N–H 92, N–CA 11, N–C 15, CA–C 55, CA–H 2, and H–C 4.
- Uses the `sphten-liouv` basis formalism with approximation `none`, then constructs the spin system with `create` and `basis`.
- Sets `J_nh=92`, `T=25e-3`, and `delta2=3e-3`. The sequence spins are `{'15N','13C','1H'}`; sweeps are `[3000 8000 2000]`, offsets are `[-6500 25000 4000]`, acquisition points are `[64 64 64]`, zero-filled sizes are `[256 256 256]`, and axis units are `ppm`.
- Simulates the sequence with `liquid(spin_system,@hncaco,parameters,'nmr')`.

## Processing and output

Applies squared-cosine apodisation in all three dimensions to the four FIDs (`pos_pos`, `pos_neg`, `neg_pos`, and `neg_neg`). Zero-filled, shifted FFTs along dimension 3 are combined as `f3_pos=f3_pos_pos-conj(f3_neg_neg)` and `f3_neg=f3_neg_pos-conj(f3_pos_neg)` to obtain the absorption part of the F3 signal. Shifted FFTs along dimension 2 are combined as `f3f2=f3f2_pos-conj(f3f2_neg)` for the absorption part of the F2 signal. A shifted FFT along dimension 1 produces `spectrum`. The imaginary part is displayed with `plot_3d` in positive mode.