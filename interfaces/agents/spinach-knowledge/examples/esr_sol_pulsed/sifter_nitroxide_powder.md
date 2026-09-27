# examples/esr_sol_pulsed/sifter_nitroxide_powder.m

- Signature: `sifter_nitroxide_powder()`

## Purpose

Powder-averaged two-dimensional SIFTER simulation for a nitroxide electron pair coupled to two ¹⁴N nuclei at 0.33 T. The example calculates the imaginary part of the time-domain signal in Liouville space and displays both the 2D SIFTER map and its diagonal. Calculation time: minutes.

## Spin system

The system has two electron spins with matching anisotropic g principal values [2.0087, 2.0058, 2.0018], and two ¹⁴N spins. The electron pair is specified at coordinates [0, 0, 0] and [0, 0, 20]. Each electron is coupled to its associated nitrogen with principal hyperfine values [19.8977, 20.1780, 102.8516] MHz; the listed tensor orientations are zero. The basis is sphten Liouville space without approximation, with `14N` longitudinal components specified.

## Simulation

The sequence starts from the electron `Lz` state, detects `L+`, and specifies `Ly` and `Lx` pulse operators. It uses zero offset, 200 points, an 8 ns timestep, and the `rep_2ang_3200pts_sph` powder grid. The resulting imaginary signal is plotted on a nanosecond time axis as a 2D image and as its diagonal trace.
