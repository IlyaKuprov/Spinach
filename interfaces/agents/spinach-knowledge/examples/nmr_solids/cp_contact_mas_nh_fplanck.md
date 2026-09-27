# examples/nmr_solids/cp_contact_mas_nh_fplanck.m

- Signature: `cp_contact_mas_nh_fplanck()`

## Purpose

Simulates a cross-polarisation contact experiment in the doubly rotating frame for a two-spin 15N-1H system under MAS, starting from the thermal-equilibrium condition. The source estimates a runtime of seconds.

## Model and calculation

- Uses a 9.394 T field, zero isotropic shifts, coordinates separated by 1.05 in the source input, and temperature 298 K. The basis is `sphten-liouv` with no approximation.
- The MAS rate is 10,000 Hz about axis `[sqrt(2/3) 0 sqrt(1/3)]`; maximum rank is 4 and the powder grid is `rep_2ang_800pts_sph`.
- Applies 100 RF steps of 10 microseconds with powers 50 kHz (1H) and 40 kHz (15N). The source assigns irradiation operators `{Hy Nx}`, excitation operators `{Hx Ny}`, and detects 15N along `Lx`.
- The source invokes `singlerot` with `@cp_contact_hard` and plots the real signal against cumulative contact time.
