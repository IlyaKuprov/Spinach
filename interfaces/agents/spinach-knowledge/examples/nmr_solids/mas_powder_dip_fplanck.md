# examples/nmr_solids/mas_powder_dip_fplanck.m

- Signature: `mas_powder_dip_fplanck()`

## Purpose

Spinning-powder pulse-acquire experiment on two dipolar-coupled protons. The source header identifies the Fokker–Planck formalism; the simulation call is `singlerot(...)`. Calculation time: seconds.

## Physical / mathematical content

- The two `1H` spins have isotropic Zeeman shifts 5.0 and -2.0 and coordinates `[0 0 0]` and `[0 3.9 0.1]`; the system is set to 14.1 T.
- The source comment cites [doi:10.1016/j.jmr.2016.07.005](https://doi.org/10.1016/j.jmr.2016.07.005). The MAS rate is 1000 Hz along `[1 1 1]`.

## Numerical / algorithmic content

- Uses the spherical-tensor Liouville-space basis with no approximation and projection +1; the angular grid is `leb_2ang_rank_17` with maximum rank 17.
- Acquires 512 points over a sweep of `2e4`, zero-fills to 4096, applies exponential apodisation parameter 6, and Fourier transforms the FID.

## Implementation structure

- Define the two-spin system and basis, build the Spinach system, set acquisition and rotor parameters, call `singlerot(spin_system,@acquire,parameters,'nmr')`, apodise, Fourier transform, and plot.
