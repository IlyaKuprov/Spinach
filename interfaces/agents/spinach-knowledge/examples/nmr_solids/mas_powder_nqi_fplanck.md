# examples/nmr_solids/mas_powder_nqi_fplanck.m

- Signature: `mas_powder_nqi_fplanck()`

## Purpose

Simulates the powder MAS spectrum of a single quadrupolar deuterium nucleus. The source header identifies Fokker–Planck theory and states that perturbative corrections to the rotating-frame transformation are not applied; the simulation call is `singlerot(...)`. Calculation time: seconds.

## Physical / mathematical content

- The source specifies one `2H` nucleus at 9.4 T, with quadrupolar coupling eigenvalues `[-1e3 -2e3 3e3]` and Euler angles `[0 0 0]`.
- The rotor axis is `[1 1 1]` at 1000 Hz, with `2H` as the initial state and detected operator.

## Numerical / algorithmic content

- Uses the spherical-tensor Liouville-space basis with no approximation and projection +1. The source sets maximum rank 17 and grid `leb_2ang_rank_17`.
- Acquires 512 points over a `2e4` sweep, zero-fills to 4096, applies exponential apodisation parameter 6, then Fourier transforms and plots the real spectrum.

## Implementation structure

- Define the quadrupolar spin system and basis, configure the experiment, call `singlerot(spin_system,@acquire,parameters,'nmr')`, apodise, Fourier transform, and plot.
