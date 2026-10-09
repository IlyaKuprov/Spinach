# examples/nmr_solids/mas_powder_csa_floquet.m

Source: [examples/nmr_solids/mas_powder_csa_floquet.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/mas_powder_csa_floquet.m)

## Model

The header calls this a powder MAS spectrum of a single anisotropically shielded proton. The active system declaration instead contains two `1H` spins and two Zeeman eigenvalue triplets, `[-2 -2 4]-5` and `[-1 -3 4]+5`, with both Euler-angle triples set to `[0 0 0]`. The model field is 14.1 T. The script specifies these shielding tensors but no dipolar coupling, quadrupolar interaction, or RF pulse sequence; it also sets an empty decoupling list.

## Calculation and display

The basis is `sphten-liouv` with no approximation and the `+1` projection. The MAS axis is `[1 1 1]` and the rotor rate is 500 Hz. The source selects the `leb_2ang_rank_17` grid and maximum rank 17, then calls `floquet(spin_system,@acquire,parameters,'nmr')`. Acquisition is on `1H` from an `L+` initial state with an `L+` receiver. The sweep is 20 kHz, with 512 acquired points, zero-filled to 4096, zero offset, ppm axis units, and inverted axis. Exponential apodisation uses parameter 6 before Fourier transformation; the plotted trace is the real spectrum. These are simulation and display settings, not an experimentally measured spectrum.
