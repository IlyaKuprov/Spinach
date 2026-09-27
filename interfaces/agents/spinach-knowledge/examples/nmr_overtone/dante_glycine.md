# examples/nmr_overtone/dante_glycine.m

- Signature: `dante_glycine()`

## Purpose

Calculates a `14N` overtone DANTE spectrum of glycine with the Fokker–Planck formalism. The glycine quadrupolar tensor is attributed to O'Dell and Ratcliffe, [Chemical Physics Letters (2011)](https://doi.org/10.1016/j.cplett.2011.08.030); the source estimates minutes of calculation time.

## Model and calculation

The single-`14N` system is at 14.1 T, with quadrupole parameters 1.18 MHz and asymmetry 0.53 (spin 1), and a scalar Zeeman entry of 32.4. The model uses diagonal damping relaxation at rate 300, the `sphten-liouv` basis without approximation, and disables Krylov and trajectory-level options.

At the magic angle, the spectrum setup uses rank 7, rate −19.840 kHz, grid `rep_2ang_1600pts_sph`, a −60 to 80 kHz sweep, and 2048 points with 2048-point zero-fill. The initial state is `14N` `Lz`; the receiver is magic-angle weighted. The average-treatment DANTE sequence has pulse amplitude 2π × 55 kHz divided by sin of the magic angle, 10 μs pulse duration, 48 kHz RF frequency, four periods, and two pulses. The code runs `singlerot` with `overtone_dante` and applies a phase factor `exp(-1i*2.12)`.
