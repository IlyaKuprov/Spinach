# examples/extremes/fluoroisooctane.m

- Signature: `fluoroisooctane()`
- Source: [`examples/extremes/fluoroisooctane.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/extremes/fluoroisooctane.m)

## Purpose and spin model

This is the source's deliberately adversarial large-spin-system example, attributed in its header to Art Bochevarov at Schrödinger Inc. The header explains the computational point: an IK-2 Liouville-space approximation would generate an exceedingly large basis, so the example instead uses Hilbert space with permutation-symmetry factorisation. The code contains 18 spins (17 `1H` and one `19F`), uses `sys.magnet=11.74` under a “Magnet induction” comment, and supplies chemical shifts and scalar couplings. The source does not annotate a unit for the magnet, shifts, or couplings.

The basis is `zeeman-hilb` with `approximation='none'` and three `S3` groups over spin sets `[1 2 3]`, `[4 5 6]`, and `[7 8 9]`. In the coupling table, the source assigns values 48.2, 23.6, and 6.5 to the central larger couplings, 6.0 to the listed neighboring proton couplings, and 1.0 to the smaller couplings it labels tert-butyl and isopropyl. These are input values as written; their units are not stated in the file.

## Proton acquisition and spectrum

The observable is a simulated `1H` NMR spectrum. The code sets both initial state and coil to proton `L+`, leaves decoupling empty, and calls `liquid(spin_system,@acquire,parameters,'nmr')`—there is no explicit pulse sequence. Acquisition settings are `offset=1290`, `sweep=2000`, `npoints=4096`, and `zerofill=8192`; the source explicitly requests a ppm axis and inverted display axis but does not label offset/sweep units. The FID is apodised with a Gaussian parameter of 6, Fourier transformed, and plotted with `plot_1d`.

## Scope and limitations

The source comment estimates a calculation time of hours and does not state a hardware configuration. This is an example of basis-size management for this model, not a measured comparison of alternative calculations; the file provides no output values, benchmark table, DOI, or experimental validation.
