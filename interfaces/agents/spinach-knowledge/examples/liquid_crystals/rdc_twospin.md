# examples/liquid_crystals/rdc_twospin.m

Source: [examples/liquid_crystals/rdc_twospin.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/liquid_crystals/rdc_twospin.m)

## System and phenomenon

A two-spin 1H–13C C–H system is used to simulate a CLIP-HSQC spectrum in a liquid crystal. The code supplies a user-defined order matrix and requests residual-dipolar-coupling support (`parameters.needs={'rdc'}`).

## Model and parameters

The spin system is set at `sys.magnet=5.9`, with chemical shifts `5.0` and `65.0`, a scalar coupling of `140`, and the order matrix `diag([1e-3 2e-3 -3e-3])`. The source does not annotate units for the field or coupling values. It uses the `sphten-liouv` formalism with no basis approximation.

The CLIP-HSQC call uses `parameters.J=140`, sweeps `[3000 1000]`, offsets `[4250 1200]`, `128` points per dimension and zero-fills to `[512 512]`; the displayed axis units are ppm. The selected spin order is `13C`, then `1H`.

## Calculation and observable

`liquid(spin_system,@clip_hsqc,parameters,'nmr')` produces positive and negative phase-cycle FIDs. Each is apodised with squared-cosine windows, Fourier transformed along F2, combined as `f1_pos+conj(f1_neg)` for States processing, then Fourier transformed along F1. The plotted observable is the real part of the resulting 2D spectrum.

## Scope

The calculation uses one fixed order matrix and one scalar coupling for this two-spin system.
