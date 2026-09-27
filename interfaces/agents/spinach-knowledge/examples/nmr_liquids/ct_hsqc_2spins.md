# examples/nmr_liquids/ct_hsqc_2spins.m

- Signature: `ct_hsqc_2spins()`

## Purpose

CT HSQC spectrum of 2 spin system. Calculation time: seconds

## Physical / mathematical content

This is a two-spin heteronuclear constant-time HSQC example with one `13C` and one `1H` site, coupled by 140.0. It uses Spinach's liquid-state `ct_hsqc` simulation and combines the positive- and negative-frequency signals using States-style quadrature reconstruction.

## Numerical / algorithmic content

The source sets field value 5.9, chemical shifts 50.00 and 3.00, and `parameters.J=140`. It uses sweeps [2500 950], offsets [3000 600], 128 points and 512 zero-fill points, with `13C` decoupling in F2. It dilutes the system into carbon isotopomers, simulates each, applies squared-cosine apodisation separately to positive and negative FIDs, Fourier-transforms in F2, combines them as `f1_pos+conj(f1_neg)`, accumulates the F1 transforms, and plots the real spectrum in negative display mode.

## Implementation structure

The function constructs the heteronuclear system, creates the carbon-diluted subsystem list and preallocates the zero-filled complex spectrum. A parallel loop builds each subsystem basis, runs `ct_hsqc`, processes both signal components and adds their F1 transforms to the accumulated spectrum.
