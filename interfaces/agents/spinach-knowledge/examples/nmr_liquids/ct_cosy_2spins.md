# examples/nmr_liquids/ct_cosy_2spins.m

- MATLAB implementation: [examples/nmr_liquids/ct_cosy_2spins.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/ct_cosy_2spins.m)

- Signature: `ct_cosy_2spins()`

## Purpose

A two-proton constant-time COSY calculation, using the chemical-shift and coupling assignment cited in [DOI 10.1002/jhet.5570250160](http://dx.doi.org/10.1002/jhet.5570250160). The source labels the run as taking minutes.

## Spin system and sequence

The model has two 1H sites at 2.00 and 5.00 ppm in a 5.9 T field, joined by a 7.0 Hz scalar coupling. The wrapper uses the complete spherical-tensor Liouville basis (`sphten-liouv`, no basis approximation) and calls the liquid-state `ct_cosy` sequence through `liquid(...,'nmr')`. It requests proton detection, sets offset to 500 and sweeps to [2000 2000], with 512 acquired points and 2048 zero-fill points per dimension; the displayed axes are requested in ppm.

The wrapper contains no separate mixing-time, phase-cycle table, or receiver-phase setting: it supplies the spin system and acquisition parameters to `ct_cosy`. Treat sequence-internal timing and coherence/receiver handling as belonging to that sequence implementation, not as independently specified here.

## Processing and output

Both dimensions receive squared-cosine apodisation. The resulting 2D FID is Fourier transformed with a shifted 2D FFT, and its magnitude is displayed with the positive plotting convention. This is a simulated spectrum from the two-site model; the example does not load experimental data or report an experimental comparison.
