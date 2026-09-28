# examples/nmr_overtone/mas_boron_2.m

- Signature: `mas_boron_2()`

## Purpose

Simulates an overtone `10B` MAS NMR spectrum with the sample spinning in the JEOL direction. The source credits parameters to Nghia Duong and Yusuke Nishiyama and specifies realistic RF power and pulse width; it estimates hours of calculation time.

## Model and calculation

The `10B` system is at 16.4 T, with quadrupole parameters 0.7 MHz, asymmetry 0, and spin 3. Diagonal damping relaxation is used at rate 100, and trajectory-level options are disabled.

At the magic angle, the sequence uses rank 12, a 70 kHz spinning rate, grid `rep_2ang_200pts_sph`, a −141 to −139 kHz sweep, and 256 points with 256-point zero-fill. The average-treatment RF settings are 2π × 50 kHz divided by sin of the magic angle, 2 ms duration, and −140 kHz frequency. The code runs `singlerot` with `overtone_pa` and applies phase factor `exp(1i*1.45)`.
