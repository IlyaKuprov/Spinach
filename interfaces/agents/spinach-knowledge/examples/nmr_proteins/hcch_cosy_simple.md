# examples/nmr_proteins/hcch_cosy_simple.m

- Signature: `hcch_cosy_simple()`

## Purpose

Simulates a 3D HCCH COSY experiment on a small protein fragment. The source notes a calculation time of seconds.

## Physical / mathematical content

- The five-spin system contains three `13C` spins (`C`, `CA`, `CB`) and two `1H` spins (`HA`, `HB`) at a magnetic field of 14.1. The scalar Zeeman values are `[180 60 40 4 2]`.
- Scalar couplings are specified for spin pairs (2,3), (1,3), (1,2), (3,5), (4,5), (3,4), (2,4), and (1,4), with values `35.0`, `0.3`, `55.0`, `130.0`, `7.0`, `6.0`, `140.0`, and `4.0`, respectively. The (5,5) entry is set to zero.
- The sequence is configured for `1H`–`13C`–`1H` dimensions, with `13C` decoupling in F3. The pulse sequence itself is supplied by `hcch_cosy` and is not defined on this page.

## Numerical / algorithmic content

- Uses the `sphten-liouv` formalism with approximation `none`. Sequence parameters are `J_ch=140`, `J_cc=35`, and `delta=1.1e-3`.
- Sweep widths are `[4000 7500 4000]`, offsets are `[2000 7000 2000]`, and acquisition sizes are `[64 64 64]`. Each dimension is zero-filled to 256 points for its Fourier transform; axis units are `ppm`.
- Applies `sqcos` apodisation in all three dimensions to each of the four simulated signal components. Fourier transforms use `fft` and `fftshift` along F3, F2, then F1. Conjugate component combinations form the absorption parts after the F3 and F2 transforms.

## Implementation structure

1. Define the spin system, interactions, basis, and sequence parameters; construct the Spinach spin system with `create` and `basis`.
2. Simulate with `liquid(spin_system,@hcch_cosy,parameters,'nmr')`, producing `pos_pos`, `pos_neg`, `neg_pos`, and `neg_neg` components.
3. Apodise all four components and transform along F3; combine them as `f3_pos=f3_pos_pos+conj(f3_neg_neg)` and `f3_neg=f3_neg_pos+conj(f3_pos_neg)`.
4. Transform both combinations along F2 and form `f3f2=f3f2_pos+conj(f3f2_neg)`; transform along F1 to obtain `spectrum`.
5. Create a figure and plot `imag(spectrum)` using `plot_3d(spin_system,imag(spectrum),parameters,10,[0.05 0.25 0.05 0.25],2,'positive')`.