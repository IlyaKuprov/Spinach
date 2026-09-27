# examples/nmr_overtone/mas_boron_3.m

- Signature: `mas_boron_3()`

## Purpose

Simulates an overtone `10B` MAS NMR spectrum with the sample spinning in the JEOL direction. The source credits parameters to Nghia Duong and Yusuke Nishiyama and describes an unphysically strong pulse to obtain a panoramic spectrum. It estimates hours of calculation time.

## Model and calculation

The `10B` system is at 16.4 T, with quadrupole parameters 0.7 MHz, asymmetry 0, and spin 3. Diagonal damping relaxation is used at rate 1000, and trajectory-level options are disabled.

The sequence setup uses rank 12, a 70 kHz spinning rate, grid `rep_2ang_200pts_oct`, and a −200 to 200 kHz sweep with 4096 points and 4096-point zero-fill. The initial state and receiver are both `10B` `Lz`; the code runs `singlerot` with `overtone_a`.
