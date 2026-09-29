# examples/nmr_liquids/ct_cosy_three_spin.m

- MATLAB implementation: [examples/nmr_liquids/ct_cosy_three_spin.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/ct_cosy_three_spin.m)

- Signature: `ct_cosy_three_spin()`

## Purpose

A compact three-spin constant-time COSY example. The source describes the calculation as taking seconds.

## Spin system and sequence

The model contains three 1H sites at 2.70, 4.10, and 6.50 ppm in a 14.1 T field. Every pair is coupled: J(1,2)=10 Hz, J(2,3)=8 Hz, and J(1,3)=4 Hz. The full `sphten-liouv` basis is used without approximation; the source also enables greedy system reduction with proximity cutoff 4.0.

The wrapper passes the system to the liquid-state `ct_cosy` sequence using `liquid(...,'nmr')`, selects 1H and sets offset to 2700 and sweeps to [3500 3500]. It requests 256 points and 512 zero-fill points in each dimension, with axis units set to ppm. It does not supply a separate mixing-time value, phase-cycle table, or receiver-phase setting, so those sequence details are not specified independently by this wrapper.

## Processing and output

The two-dimensional FID is squared-cosine apodised along both dimensions, transformed with a shifted 2D FFT, and plotted as magnitude in positive mode. This example demonstrates the response of the specified three-site coupling network; the source contains no experimental acquisition or experimental comparison.
