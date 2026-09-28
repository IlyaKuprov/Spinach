# examples/esr_sol_pulsed/holeburn_gd_dota_powder.m

- Signature: `holeburn_gd_dota_powder()`

## Purpose

Simulates Gd(III) spectral hole burning with a soft pulse. It samples a zero-field-splitting (ZFS) distribution, performs powder averaging, and uses a numerical second-order rotating-frame transformation.

## Model and calculation

For each weighted ZFS sample returned by `zfs_sampling(30,5,1e-2)`, the script builds an `E8` system at 3.5 T in a spherical-tensor Liouville basis. It computes spectra with the pulse power set first to zero and then to `2π × 10⁷ rad s⁻¹`, using a 50 ns, rank-2 pulse at −0.5 GHz. Both signals are apodised and Fourier transformed, then accumulated with the sample weights for plotting. The source notes that non-central Gd(III) transition holes are very shallow; estimated run time is minutes.

The ZFS-distribution parameters are taken from Fig. 5 of Raitsimring et al., *Applied Magnetic Resonance* **28**, 281–295 (2005).
