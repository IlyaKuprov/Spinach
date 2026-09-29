# examples/nmr_solids/mas_powder_csa_gridfree.m

Source: [examples/nmr_solids/mas_powder_csa_gridfree.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/mas_powder_csa_gridfree.m)

## Model

This source models two `1H` spins at 14.1 T with Zeeman eigenvalue triplets `[-2 -2 4]-5` and `[-1 -3 4]+5`, each with Euler angles `[0 0 0]`. It specifies anisotropic shielding tensors, but no dipolar coupling, quadrupolar interaction, or RF pulse sequence. The header describes a powder MAS spectrum using a grid-free Fokker-Planck formalism.

## Calculation and display

The basis is `sphten-liouv` with no approximation and the `+1` projection. The rotor axis is `[1 1 1]` at 500 Hz. The experiment settings acquire on `1H` from and to `L+` states, with an empty decoupling list, a 20 kHz sweep, 512 points, and zero-fill to 4096. No explicit orientation grid is set. The active call is `gridfree(spin_system,@acquire,parameters,'nmr')`.

After exponential apodisation with parameter 6, the code Fourier transforms the FID and plots the real spectrum. The axis units are ppm and the axis is inverted. These are simulation and display settings, not an experimentally measured spectrum.
