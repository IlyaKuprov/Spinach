# examples/nmr_proteins/hcanh_simple.m

- Signature: `hcanh_simple()`

## Purpose

A minimal example of H(CA)NH pulse sequence simulation. Calculation time: seconds.

## Physical / mathematical content

The example defines a five-spin liquid-state NMR system at a magnetic field of 14.1, with isotopes `{'15N','13C','1H','13C','1H'}` labelled `{'N','CA','H','C','HA'}`. Scalar Zeeman values are `{110 60 8 180 4}`. The specified scalar couplings, indexed by spin number, are (1,3) = 92, (1,2) = 11, (1,4) = 15, (1,5) = 1, (2,4) = 55, (2,3) = 2, (2,5) = 140, (3,4) = 4, (3,5) = 8, and (4,5) = 4.

The calculation uses the `sphten-liouv` basis formalism with `approximation='none'`. It simulates the `@hcanh` sequence through `liquid(...,'nmr')` with dimensions `{'1H','15N','1H'}`, sweep widths `[5000 3000 5000]`, offsets `[3600 -6600 3600]`, `npoints=[64 64 64]`, `zerofill=[256 256 256]`, and axis units `ppm`.

## Numerical / algorithmic content

The four returned coherence components (`pos_pos`, `pos_neg`, `neg_pos`, and `neg_neg`) each receive squared-cosine (`sqcos`) apodisation in all three dimensions. Shifted, zero-filled FFTs are applied along F3, then F2, then F1. The F3 components are combined as `f3_pos=f3_pos_pos+conj(f3_neg_neg)` and `f3_neg=f3_neg_pos+conj(f3_pos_neg)`; after the F2 transforms, `f3f2=f3f2_pos+conj(f3f2_neg)`. The final F1 transform produces `spectrum`.

## Implementation structure

1. Define the magnetic field, spin isotopes and labels, Zeeman values, and scalar couplings.
2. Create the spin system and construct its basis with `create` and `basis`.
3. Set the three-dimensional sequence and acquisition parameters, then run `liquid(spin_system,@hcanh,parameters,'nmr')`.
4. Apodise the four FID components; perform the F3, F2, and F1 FFTs and combine the specified components for absorption-mode processing.
5. Create a figure and plot `real(spectrum)` using `plot_3d(spin_system,real(spectrum),parameters,10,[0.05 0.5 0.05 0.5],2,'positive')`.