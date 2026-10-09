# examples/nmr_solids/cp_contact_mas_nh_floquet.m

Source: [examples/nmr_solids/cp_contact_mas_nh_floquet.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/cp_contact_mas_nh_floquet.m)

## Purpose and spin model

This is a simulated 15N-1H powder CP contact experiment in the doubly rotating frame, starting from the isotropic-equilibrium condition. The two-spin model uses `sys.magnet=9.394`, zero isotropic shifts for both spins, coordinates `[0 0 0]` and `[0 0 1.05]`, and temperature value 298. The file does not state a unit for the coordinate separation or temperature. It uses the full `sphten-liouv` basis without approximation. No experimental FID or spectrum is read as input.

## MAS contact sequence and plotted signal

The MAS settings are rate 10000, axis `[sqrt(2/3) 0 sqrt(1/3)]`, maximum rank 4, and powder grid `rep_2ang_800pts_sph`. The rate value has no unit attached in the source. The spin order is `{'1H','15N'}`; across 100 steps the irradiation powers are 5e4 for 1H and 4e4 for 15N, with irradiation operators `{Hy Nx}` and excitation operators `{Hx Ny}`. Each time step is 1e-5, and the source labels the plotted cumulative time axis in seconds. It requests `iso_eq` and detects on the 15N `Lx` state.

The sequence is propagated by `floquet` with `@cp_contact_hard` in NMR mode. The program plots the real signal against cumulative contact time and labels the ordinate as the 15N `S_X` expectation value. This is a contact-time transfer trace, not a frequency-domain spectrum; the file supplies a simulated curve rather than measured data. The source estimates a runtime of seconds and does not report numerical transfer values or an experimental comparison.
