# examples/liquid_crystals/rdc_fourspin.m

Source: https://github.com/IlyaKuprov/Spinach/blob/main/examples/liquid_crystals/rdc_fourspin.m

- Signature: `rdc_fourspin()`

## Physical model

This seconds-scale example simulates a CLIP-HSQC spectrum for a four-spin system containing two 1H and two 13C spins at 5.9 T. Scalar-shift inputs are `[1, 10, 80, 50]`; the nonzero scalar couplings are `J12=2.5`, `J13=250`, `J14=50`, and `J23=50`. The four coordinates lie along the x axis at `[1, 2, 3, 4]` in the source's coordinate units. A user-specified traceless order matrix, `diag([1e-3, 2e-3, -3e-3])`, provides the anisotropic ordering that makes residual dipolar coupling relevant. The basis is sphten-liouv with no approximation.

## Sequence and spectrum

The sequence call is `liquid(spin_system,@clip_hsqc,parameters,'nmr')`, with `parameters.needs={'rdc'}`. It sets both sweep widths to 5000 Hz, offsets to `[4250, 1200] Hz`, 128 points per dimension, zero filling to 512 per dimension, spins ordered as 13C then 1H, and `J=140` Hz. The positive and negative FIDs receive square-cosine apodisation in both dimensions. The script Fourier transforms the first dimension, forms the States signal as `f1_pos+conj(f1_neg)`, transforms the second dimension, and plots the real spectrum in ppm.

The order matrix and four-spin parameters are explicit inputs to this demonstration; the source does not claim a fitted molecular structure or a particular experimental spectrum.
