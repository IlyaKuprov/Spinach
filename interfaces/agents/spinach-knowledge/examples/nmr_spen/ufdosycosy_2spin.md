# examples/nmr_spen/ufdosycosy_2spin.m

- Signature: `ufdosycosy_2spin()`

## Purpose

Simulates ultrafast 3D DOSY-COSY for a two-spin system. The source estimates hours on an NVIDIA Tesla A100 and much longer on CPU. Authors: Ludmilla Guduff and Jean-Nicolas Dumez.

## Model and sequence

The system contains two 1H spins at 14.1 T, with chemical shifts 2.0 and -1.3 and a 15 Hz scalar coupling. The basis is sphten-liouv without approximation. Redfield relaxation phantoms are empty; flow is zero and the diffusion coefficient is 4e-10 m^2/s. The sample is 0.015 m long with 3000 points and a period-7 derivative.

The imaging call uses spendosycosy. Acquisition is 0.100 s with 1.5e-6 s dwell time, 128 and 64 points in the two listed dimensions, 256 loops, sweep 1302 Hz, gradient Ga=0.52 T/m, and offset 600 Hz. Encoding uses 1000 pulse points, smfactor 0.1, Te=0.0015 s, Tau=0.0016 s, BW=110000 Hz, Ge=0.2535 T/m, and smoothed chirps. The fit settings are dscale=1e-10, fovmin=-0.0041, fovmax=+0.0041, Hamming apodization, a sine window, and the keeler_corr model.

## Processing

The returned signal is Fourier transformed along three dimensions and displayed as the magnitude raised to the one-half power.
