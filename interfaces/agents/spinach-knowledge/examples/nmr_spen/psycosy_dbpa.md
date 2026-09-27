# examples/nmr_spen/psycosy_dbpa.m

- Signature: `psycosy_dbpa()`

## Purpose

PSYCOSY of DBPA (dibromopropionic acid) ring. Calculation time: minutes on NVidia Tesla A100, much longer on CPU

## Physical / mathematical content

- The source defines the DBPA spin system through its isotope, chemical-shift, and scalar-coupling data.
- It runs the `@psycosy` sequence using Spinach's imaging function; the sequence settings are defined in the source.

## Numerical / algorithmic content

- The simulated 2D FID is square-sine apodised along both dimensions and transformed by a 2D FFT. The script plots the magnitude of that spectrum.

## Implementation structure

- Creates the spin system and basis, sets sequence parameters, runs imaging, and processes and plots the spectrum.
