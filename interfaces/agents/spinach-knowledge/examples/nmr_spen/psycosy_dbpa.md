# examples/nmr_spen/psycosy_dbpa.m

- Signature: `psycosy_dbpa()`

## Purpose

PSYCOSY of DBPA (dibromopropionic acid) ring. Calculation time: minutes on NVidia Tesla A100, much longer on CPU

## Physical / mathematical content

- The source defines a four-1H DBPA spin system at 14.1 T, with chemical shifts 4.49, 3.9, 3.7, and 4.2 ppm and nonzero scalar couplings of 11.3, 10.1, and 4.3 Hz. It models a 15 mm sample with 100 spatial points.
- It runs `@psycosy` through Spinach's imaging function. The 2D acquisition uses 512 points per dimension, zero-filled to 1024, with a 600 Hz sweep and 2460 Hz offset; sequence settings include a 25 ms mixing time, 0.01 T/m gradient, and a 20° saltire chirp (15 ms pulse, 10 kHz sweep, 50 ms gradient duration).

## Numerical / algorithmic content

- The simulated 2D FID is square-sine apodised along both dimensions and transformed by a 2D FFT. The script plots the magnitude of that spectrum.

## Implementation structure

- Creates the spin system and basis, sets sequence parameters, runs imaging, and processes and plots the spectrum.
