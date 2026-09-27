# examples/nmr_solids/mas_powder_dip_floquet.m

- Signature: `mas_powder_dip_floquet()`

## Purpose

Simulates a two-proton spinning-powder pulse-acquire experiment with dipolar coupling using Floquet theory. Calculation time: seconds.

## Physical / mathematical content

- The two `1H` spins have isotropic Zeeman shifts 5.0 and -2.0 and coordinates `[0 0 0]` and `[0 3.9 0.1]`; the system is set to 14.1 T.
- The rotor axis is `[1 1 1]` at 1000 Hz. The simulation obtains the acquisition FID with `floquet(spin_system,@acquire,parameters,'nmr')`.

## Numerical / algorithmic content

- Uses the spherical-tensor Liouville-space basis with no approximation and projection +1. The MAS grid is `leb_2ang_rank_17`, with maximum rank 17.
- Acquires 512 points over a sweep of `2e4`, zero-fills to 4096, applies exponential apodisation parameter 6, then Fourier transforms and plots the real spectrum.

## Implementation structure

- Define the two-spin system and basis, build the Spinach system, set pulse-acquire and MAS parameters, run the Floquet acquisition, apodise, Fourier transform, and plot.
