# examples/nmr_solids/mas_powder_nqi_gridfree.m

- Signature: `mas_powder_nqi_gridfree()`

## Purpose

Powder magic angle spinning spectrum of a single quadrupolar deuterium nucleus using grid-free Fokker-Planck MAS formalism. Second order corrections to the rotating frame transformation are not applied. Calculation time: minutes

## Physical / mathematical content

- Models a single `2H` nucleus with a quadrupolar coupling tensor whose eigenvalues are `[-1e3 -2e3 3e3]` and Euler angles `[0 0 0]`.
- Simulates powder MAS with the grid-free Fokker-Planck formalism; the example explicitly states that second-order rotating-frame transformation corrections are not applied.

## Numerical / algorithmic content

- Acquires the FID with `gridfree`, applies exponential apodisation with parameter 6, zero-fills to 4096 points, and plots the real Fourier-transformed spectrum.

## Implementation structure

- Defines a single-deuterium system at 9.4 T with quadrupolar coupling eigenvalues `[-1e3 -2e3 3e3]` and zero Euler angles.
- Uses the `sphten-liouv` basis without approximation and retains the `+1` projection.
- Sets the MAS axis to `[1 1 1]`, rate to `1e3`, maximum rank to 17, sweep to `2e4`, and acquisition to 512 points with 4096-point zero filling.
- Creates the initial state and receiver as `L+` on `2H`, then runs `gridfree` with `acquire` in NMR mode.
- Applies exponential apodisation with parameter 6, Fourier transforms the FID, and plots the real spectrum.
