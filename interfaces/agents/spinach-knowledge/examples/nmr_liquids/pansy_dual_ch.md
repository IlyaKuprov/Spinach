# examples/nmr_liquids/pansy_dual_ch.m

- Signature: `pansy_dual_ch()`
- Source: [examples/nmr_liquids/pansy_dual_ch.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/pansy_dual_ch.m)

## What it models and loads

A two-channel PANSY-COSY simulation for camphor with natural 13C content, as identified by the source comment. The wrapper reads `../standard_systems/camphor.log`, described as vacuum-DFT input; its comment attributes coordinates, shieldings and scalar couplings to DFT. This is a computed-chemistry input, not measured NMR data. The source's estimate of seconds is a comment, not a measured run time.

## Spin system

The DFT data are passed to `g2spinach` with the isotope mapping H to 1H and C to 13C, positional values [31.5, 189.2], and options `min_j=3.0` and `no_xyz=1`. The wrapper does not name or unit-label the two positional values. It sets the field to 14.1 (unit not stated) and uses a spherical-tensor Liouville basis, IK-2 approximation, scalar-coupling connectivity and proximity level 1. `dilute(spin_system,'13C')` supplies the set of isotopomer subsystems; the wrapper simulates each in a parallel loop and accumulates the resulting channels. It does not state an experimental abundance measurement.

## Acquisition and processing

The sequence callback is `pansy_cosy`; the wrapper sets sweeps [1800, 9000], offsets [900, 4500], 128 points per dimension, zero filling [512, 512], spins [1H, 13C], and ppm axis units. It applies squared-cosine apodisation in both dimensions to each of the two returned signal components, `fid.aa` and `fid.ab`, then performs 2D Fourier transforms and plots their magnitudes. Sweep and offset units are not stated. Pulse, gradient, receiver and channel-generation internals are delegated to `pansy_cosy` and are not inferred from this wrapper.

## Output and limits

The wrapper plots a homonuclear-side spectrum from the accumulated `aa` component and a heteronuclear-side spectrum from `ab`; it does not write spectra to files or report measured observables. It does not specify normalisation of the accumulated isotopomer spectra.
