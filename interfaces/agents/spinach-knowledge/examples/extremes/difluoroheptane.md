# examples/extremes/difluoroheptane.m

- Signature: `difluoroheptane()`

## Purpose

19F NMR spectrum of anti-3,4-difluoroheptane (16 spins) by explicit time-domain evolution in Liouville space. WARNING: needs 32 CPU cores, 128 GB of RAM and a Titan V or later. Run time on the above: minutes

## Physical / mathematical content

- The 16-spin anti-3,4-difluoroheptane system is observed by 19F NMR at 11.7464 T; proton–fluorine couplings contribute to the fluorine spectrum.
- The computed observable is a time-domain 19F free-induction decay, converted to a frequency-domain spectrum.

## Numerical / algorithmic content

- The script uses the `sphten-liouv` formalism with greedy parallelisation, evolves the FID in time, applies exponential apodisation, and Fourier transforms it to the spectrum.

## Implementation structure

- 19F NMR spectrum of anti-3,4-difluoroheptane (16 spins) by
- explicit time-domain evolution in Liouville space.
- WARNING: needs 32 CPU cores, 128 GB of RAM and
- a Titan V or later.
- Run time on the above: minutes
- Magnet induction
- Isotopes
- Shifts
- J-couplings
- Basis set
- Greedy parallelisation
- Spinach housekeeping
