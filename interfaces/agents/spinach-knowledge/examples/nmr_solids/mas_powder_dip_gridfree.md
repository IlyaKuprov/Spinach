# examples/nmr_solids/mas_powder_dip_gridfree.m

- Signature: `mas_powder_dip_gridfree()`

## Purpose

Computes a powder MAS spectrum for two dipole-coupled protons with the grid-free Fokker–Planck method. Calculation time: minutes.

## Physical / mathematical content

- The system has two `1H` spins, isotropic Zeeman shifts 5.0 and -2.0, coordinates `[0 0 0]` and `[0 3.9 0.1]`, and a 14.1 T field.
- The rotor axis is `[1 1 1]` at 1000 Hz. The initial state and detected operator are both `L+` on `1H`.

## Numerical / algorithmic content

- Uses the spherical-tensor Liouville-space basis with no approximation and projection +1; maximum rank is 15.
- The grid-free acquisition returns 512 points over a `2e4` sweep, zero-filled to 4096. The code applies exponential apodisation parameter 6 before Fourier transformation.

## Implementation structure

- Define the two-spin system and basis, build the Spinach system, configure acquisition and MAS parameters, run `gridfree(spin_system,@acquire,parameters,'nmr')`, apodise, Fourier transform, and plot.
