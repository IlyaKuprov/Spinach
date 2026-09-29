# examples/nmr_liquids/ct_hsqc_2spins.m

- MATLAB implementation: [examples/nmr_liquids/ct_hsqc_2spins.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/ct_hsqc_2spins.m)

- Signature: `ct_hsqc_2spins()`

## Purpose

A two-spin heteronuclear constant-time HSQC simulation. One 13C site and one 1H site are coupled directly, providing a small model for following heteronuclear transfer and quadrature-sensitive processing. The source estimates a calculation time of seconds.

## Spin system and sequence

The system is at 5.9 T, with source chemical shifts of 50.00 ppm for 13C and 3.00 ppm for 1H and a 140.0 Hz scalar coupling. The sequence parameter `J` is also set to 140. The wrapper calls `liquid(...,@ct_hsqc,...,'nmr')`, requests spins [13C,1H], and sets `decouple_f2={'13C'}`. Its acquisition parameters are sweep [2500 950], offset [3000 600], 128 points and 512 zero-fill points per dimension, with ppm display axes.

The example calls `dilute(spin_system,'13C')` and processes the resulting subsystem list in a parallel loop. Both `fid.pos` and `fid.neg` are squared-cosine apodised. Their separately transformed signals are combined as `f1_pos + conj(f1_neg)` (the source's States reconstruction), then transformed along the other dimension. Thus the quadrature combination is explicit; this wrapper does not specify a separate phase-cycle table or mixing-time parameter.

## Processing and output

The accumulated spectrum is plotted using its real part and the negative plotting convention. The axes are requested in ppm. This is a calculated spectrum for the specified two-site Hamiltonian and processing choices; the script does not load or compare an experimental HSQC data set.
