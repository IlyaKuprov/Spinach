# examples/nmr_liquids/ecosy_strychnine.m

- Signature: `ecosy_strychnine()`

## Purpose

E.COSY spectrum of strychnine. Calculation time: minutes

## Physical / mathematical content

- Two-dimensional 1H E.COSY simulation of strychnine, using scalar-coupling evolution to produce correlated cross peaks.
- The cosine and sine FIDs are squared-cosine apodised. After the F2 transforms, the States-like signal is assembled as `real(f1_cos)+1i*imag(f1_sin)` and Fourier transformed along F1.

## Numerical / algorithmic content

- Uses the sphten-liouv / IK-2 basis with scalar-coupling connectivity and proximity level 1, and greedy settings with `prox_cutoff=4.0`. The field is 5.9 T; 1H acquisition settings are offset `1200`, sweep `2200`, `npoints=[512 512]`, and `zerofill=[2048 2048]`. The implementation simulates the full spin system directly, without an isotopomer loop.

## Implementation structure

- Create the 1H strychnine spin system at 5.9 T and construct the selected basis.
- Simulate E.COSY with the stated acquisition settings, apply squared-cosine windows, form the signal from the real cosine and imaginary sine components after F2 transformation, Fourier transform F1, and plot the real spectrum.
