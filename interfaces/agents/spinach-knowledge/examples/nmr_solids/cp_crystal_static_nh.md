# examples/nmr_solids/cp_crystal_static_nh.m

- Signature: cp_crystal_static_nh()

## Purpose

A simulated static single-crystal 1H-15N cross-polarisation (CP) contact curve in the doubly rotating frame. The source estimates a calculation time of seconds. This is a single specified crystal orientation, not a powder average.

## Spin model and experiment

The model contains one 15N and one 1H, with zero isotropic Zeeman shifts, coordinates [0, 0, 0] and [0, 0, 1.05] (the existing page identifies their separation as 1.05 angstrom), and temperature set to 298. The basis is sphten-liouv with no approximation. The example lists shifts and coordinates rather than explicit coupling tensors.

There is no rotor or powder grid. The selected crystal orientation is [pi/3, pi/4, pi/5]. The example requests aniso_eq for 15N, detects the 15N Lx state, and applies 5e4 Hz (50 kHz) spin-lock nutation frequency on each channel. It propagates 100 steps of 1e-5 seconds each. The specified excitation operators are Hx on 1H and Ly on 15N; the spin-lock operators are Hy on 1H and Lx on 15N.

The simulation calls crystal with the generic cp_contact_hard experiment function. That function uses the supplied ideal pi/2 excitation, spin-lock operators and powers, and time steps to compute the contact curve; it does not define the spin model or orientation. This example is CP, not HMQC, and has no DOR rotor scheme.

## Output and interpretation

The returned simulated signal is plotted as its real part against cumulative time in seconds, with the axis labelled as the 15N SX expectation value. It is not a measured spectrum or a validation result.

## Source

https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/cp_crystal_static_nh.m
https://github.com/IlyaKuprov/Spinach/blob/main/experiments/cp_contact_hard.m
