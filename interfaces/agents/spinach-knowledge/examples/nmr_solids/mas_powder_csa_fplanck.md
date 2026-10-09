# examples/nmr_solids/mas_powder_csa_fplanck.m

Source: [examples/nmr_solids/mas_powder_csa_fplanck.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/mas_powder_csa_fplanck.m)

## Model

The header describes a powder MAS spectrum of two anisotropically shielded protons using a Fokker-Planck formalism. The active code declares two `1H` spins at 14.1 T with Zeeman eigenvalue triplets `[-2 -2 4]-5` and `[-1 -3 4]+5`; both Euler-angle triples are `[0 0 0]`. No dipolar coupling, quadrupolar interaction, or RF pulse sequence is specified.

## Calculation and display

The basis is `sphten-liouv` with no approximation and the `+1` projection. The rotor axis is `[1 1 1]` at 500 Hz. Maximum rank is 17 and the orientation grid is `leb_2ang_rank_17`. Importantly, the function call in this file is `singlerot(spin_system,@acquire,parameters,'nmr')`, not a call named `fplanck` or `gridfree`; the header's formalism label and the active entry point are distinct source facts.

The signal is acquired on `1H` from and to `L+` states. The sweep is 20 kHz, with 512 points and zero-fill to 4096. The source applies exponential apodisation with parameter 6, Fourier transforms the FID, and plots the real spectrum. It does not set offset, axis units, or axis inversion in this file. The settings define a simulation and processing pipeline, not an experimentally measured spectrum.
