# examples/nmr_liquids/noesy_sucrose.m

- Signature: `noesy_sucrose()`

## Purpose

Simulates a two-dimensional proton NOESY spectrum of sucrose. The proton spin model is parsed from the vacuum-DFT sucrose log; `options.min_j=1.0` is passed to the parser as the minimum retained coupling threshold. The example does not specify units for that threshold. Its header estimates minutes for calculation time.

## Spin model and sequence

The field setting is `5.9`. The basis uses the `sphten-liouv` formalism, `IK-2` approximation, scalar-coupling connectivity, and proximity level `3`. Redfield relaxation uses IME equilibrium, temperature `298`, `rlx_keep='kite'`, and `tau_c={200e-12}`. The script enables the greedy algorithm, disables Krylov, and sets the proximity cutoff to `4.0`.

The NOESY mixing time is `0.5`; offset is `800`; sweep is `[1700 1700]`; acquired points are `[512 512]`; and zero-fill sizes are `[2048 2048]`. The selected spins are `{'1H'}`, axes are in ppm, and the simulation requests equilibrium density with `needs={'rho_eq'}`. The example leaves units unstated for mixing time, offset, sweep, and minimum-coupling threshold.

## Propagation and processing

The script calls `liquid(...,@noesy,...,'nmr')`. It applies `sqcos` apodisation to the cosine and sine FIDs in both dimensions, combines them as the States signal, Fourier-transforms F2 and F1 with the configured zero-fill lengths, and plots the negative real spectrum. The source supplies no numerical peak intensities or cross-relaxation rates.

## Source

[MATLAB example](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/noesy_sucrose.m)
