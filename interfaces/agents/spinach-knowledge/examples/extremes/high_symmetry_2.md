# examples/extremes/high_symmetry_2.m

- MATLAB implementation: [examples/extremes/high_symmetry_2.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/extremes/high_symmetry_2.m)

- Entry point: `high_symmetry_2()`.
- Source: [`examples/extremes/high_symmetry_2.m`](../../../../../examples/extremes/high_symmetry_2.m).

## Intent and spin system

This example computes a one-dimensional 31P NMR spectrum for the large, highly symmetric system described in the source comments as containing two tert-butyl groups; the system was supplied by Eberhard Matern. The isotope list defines 22 spins: two 31P and 20 1H. The script assigns the magnetic-induction parameter as 9.39798 and specifies the Zeeman and scalar-coupling tables directly. The source does not annotate units for those table entries or the field value.

The basis is `zeeman-hilb` with `approximation='none'`. The source describes the calculation as brute-force time propagation in Hilbert space; it is therefore a deliberately direct calculation of the full specified spin system, not a reduced-basis example.

## Acquisition and spectrum

The selected spin is `31P`; the initial state and receiver are both `L+`, and `decouple={}`. There is no explicit RF-pulse schedule in this script: it calls `liquid(spin_system,@acquire,parameters,'nmr')`. The scan parameters are offset `-7100`, sweep `1000`, `4096` points and zero-fill to `32768`; these numeric offset/sweep values are not unit-labelled in the source, while the displayed axis is set to `ppm` and inverted. Exponential apodisation uses parameter `5`. The code Fourier-transforms the FID with the zero-fill length, takes its real part, and plots one 1D spectrum.

## Practical limit

The source warns that the calculation needs 32+ CPU cores and 128+ GB RAM and takes hours on such a machine.
