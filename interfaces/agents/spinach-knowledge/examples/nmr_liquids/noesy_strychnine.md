# examples/nmr_liquids/noesy_strychnine.m

- Signature: `noesy_strychnine()`

## Purpose

Simulates a two-dimensional proton NOESY spectrum of strychnine. The system is supplied by `strychnine({'1H'})`, so this example selects the 1H channel. Its header estimates minutes for calculation time.

## Spin model and sequence

The field setting is `5.9`. The basis uses the `sphten-liouv` formalism, `IK-2` approximation, scalar-coupling connectivity, and proximity level `3`. Redfield relaxation uses IME equilibrium, temperature `298`, `rlx_keep='kite'`, and `tau_c={200e-12}`. The script enables the greedy algorithm, disables Krylov, and sets the proximity cutoff to `4.0`.

The NOESY mixing time is `0.5`; offset is `1200`; sweep is `[2500 2500]`; acquired points are `[512 512]`; and zero-fill sizes are `[2048 2048]`. The selected spins are `{'1H'}`, axes are in ppm, and the simulation requests the equilibrium density with `needs={'rho_eq'}`. The code does not label units for the mixing-time, offset, or sweep values.

## Propagation and processing

The script calls `liquid(...,@noesy,...,'nmr')`. It applies `sqcos` apodisation to the cosine and sine FIDs in both dimensions, combines them as the States signal, Fourier-transforms F2 and F1 with the configured zero-fill lengths, and plots the negative real spectrum. The source configures a simulated spectrum but supplies no numerical peak intensities, cross-relaxation rates, or peak-sign assignments.

## Source

[MATLAB example](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/noesy_strychnine.m)
