# examples/nmr_overtone/mas_boron_1.m

- Signature: `mas_boron_1()`

## Purpose

Simulates an overtone Z-detection `10B` MAS NMR spectrum, with the sample spinning in the JEOL direction. The source credits parameters to Nghia Duong and Yusuke Nishiyama and focuses on the most intense of the five overtone spinning sidebands. It estimates hours of calculation time.

## Model and calculation

The `10B` system is at 16.4 T, with quadrupole parameters 0.7 MHz, asymmetry 0, and spin 3. Diagonal damping relaxation is used at rate 50, and the code disables trajectory-level options.

The sequence setup uses rank 12, a 70 kHz spinning rate, grid `rep_2ang_800pts_sph`, and a −141 to −139 kHz sweep with 256 points and 256-point zero-fill. The initial state and receiver are both `10B` `Lz`; the simulation uses `singlerot` with `overtone_a`.
