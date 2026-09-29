# examples/nmr_liquids/clip_hsqc_camphor.m

- MATLAB implementation: [examples/nmr_liquids/clip_hsqc_camphor.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/clip_hsqc_camphor.m)

- Signature: `clip_hsqc_camphor()`

## Purpose and molecular spin system

This example simulates and plots a natural-abundance ¹³C CLIP-HSQC spectrum of camphor. The source estimates minutes of calculation time. It parses the vacuum DFT output at ../standard_systems/camphor.log with hydrogen and carbon-13 particle mappings and passes reference shieldings [31.8, 182.1] to `g2spinach` for the H and C channels. Its options are `min_j=3.0` and `no_xyz=0`. The example comments attribute coordinates, shielding anisotropies, and scalar couplings to DFT, then replace isotropic shifts for spins 1:26 with experimental values (ppm): 29.70, 26.80, 44.20, 59.70, 42.80, 218.1, 49.70, 17.80, 17.30, 7.80, 1.35, 1.67, 1.97, 1.31, 1.96, 0.85, 0.85, 0.85, 0.99, 0.99, 0.99, 0.90, 0.90, 0.90, 2.33, and 1.76.

The simulation uses the spherical-tensor Liouville formalism, IK-2 basis, scalar-coupling connectivity, proximity level 1, greedy algorithm setting, and proximity cutoff 4.0. The field setting is 14.1. The spin system is expanded into ¹³C isotopomers with `dilute(spin_system,'13C')` to represent natural carbon-13 abundance.

## CLIP-HSQC sequence and coherence handling

The example calls `liquid(subsystem,@clip_hsqc,parameters,'nmr')` for each isotopomer. The sequence uses F1 = ¹³C and F2 = ¹H, with active scalar coupling 140 Hz, sweep widths [8000, 1500] Hz, offsets [4000, 1000] Hz, 128 × 128 acquired points, and 512 × 512 zero filling; displayed axis units are ppm. In the sequence implementation, the initial state is proton `Lz`, detection is proton `L+`, and the coupling delay is `abs(1/(2*J))` (about 3.57 ms at J = 140 Hz). The pulse train transfers coherence through the heteronuclear coupling; the F1 evolution uses midpoint refocusing, and coherence filters select the carbon ±1 pathways with proton coherence order zero. The sequence forms the positive and negative States quadrature FIDs by pairing forward density-operator evolution with backward-propagated receiver states. The pulse-sequence source cites the CLIP-HSQC paper at https://doi.org/10.1016/j.jmr.2008.03.009.

## Processing and plotted result

The positive and negative FIDs each receive cosine-squared apodisation in both dimensions. The code Fourier transforms the direct dimension, combines the components as a States signal, Fourier transforms the indirect dimension, sums the isotopomer spectra, and plots the real part with `plot_2d`. The output is a calculated contour plot, not an experimental spectrum or a reported fit. The example provides no line-shape comparison, error metric, or claim that the simulated intensities reproduce a measured result.
