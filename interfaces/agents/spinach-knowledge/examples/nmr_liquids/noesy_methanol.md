# examples/nmr_liquids/noesy_methanol.m

- Signature: `noesy_methanol()`

## Purpose

Simulates a two-dimensional proton NOESY spectrum for the methanol spin system built from the vacuum-DFT log. The isotope map retains 13C and 1H, then the script removes the terminal OH proton, leaving the three methyl protons and 13C. The source identifies the J-couplings as values from Pecul and Helgaker and the CSA tensors as DFT-derived; it specifically notes cross-peaks between the two 13C-coupled doublet components. The field setting is `14.1`. The example header estimates seconds for calculation time.

## Spin model and sequence

All four chemical shifts are placed on resonance. The script assigns scalar-coupling values of `141` between 13C and each of the three protons, and `-11` between each proton pair; the example does not label units for these assignments. Redfield relaxation is used with IME equilibrium, temperature `298`, `rlx_keep='kite'`, and `tau_c={50e-12}`. The basis is the full `sphten-liouv` formalism with `approximation='none'`. The NOESY simulation requests equilibrium density through `needs={'rho_eq'}` and uses `parameters.spins={'1H'}`.

The sequence settings are mixing time `0.5`, offset `0`, sweep `[300 300]`, `[256 256]` acquired points, and `[1024 1024]` zero-fill sizes. The axis units are ppm. These values are reproduced as source literals; the example does not attach units to the mixing-time, offset, or sweep assignments.

## Propagation and processing

`liquid(...,@noesy,...,'nmr')` propagates the NOESY sequence. The cosine and sine FIDs are each apodised with `sqcos` in both dimensions. The script forms the States signal as cosine minus `i` times sine, Fourier-transforms the indirect dimension and then the direct dimension using the specified zero-fill lengths, and plots the negative real part. The source specifically mentions cross-peaks between the 13C-coupled doublet components but supplies no numerical cross-peak intensities or cross-relaxation rates.

## Source

[MATLAB example](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/noesy_methanol.m)
