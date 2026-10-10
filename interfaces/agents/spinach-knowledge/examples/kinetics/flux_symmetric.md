# examples/kinetics/flux_symmetric.m

- MATLAB implementation: [examples/kinetics/flux_symmetric.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/kinetics/flux_symmetric.m)

- Callable as the no-argument MATLAB function `flux_symmetric()`; it constructs the Spinach system and basis before acquisition.

## Purpose and model

This is a two-site intermolecular magnetisation-exchange simulation for two `1H` environments. The spins are separate substances undergoing additive `A+B -> A+B` replacement with swapped atom matching. Concentrations `[1 1]` and event rate `2e3` reproduce the former directed rates: 2000 in both directions. The initial `L+` state is weighted once by these concentrations; detection is unweighted. Both concentrations are invariant, so the acquisition callback freezes the additive generator at `unit_state` before ordinary linear propagation. No cross-molecular correlations are represented.

## Acquisition and observable

The full `sphten-liouv` basis is used (`bas.approximation=none`). The function passes the system, frozen-replacement acquisition callback, and NMR parameters to `liquid`, with the `1H` `L+` coil and no decoupled spins. Acquisition uses `offset=900`, `sweep=5000`, 512 points, and zero filling to 1024; its plotted axis is labelled in ppm and inverted. The FID is exponentially apodised with parameter 6, transformed with a shifted FFT, and its real spectrum is plotted with `plot_1d`.

## Scope

The source header estimates a calculation time of seconds; this is not a measured runtime. The file specifies a symmetric flux setup rather than a fitted exchange result, and it contains no numerical spectrum from which to report peak positions or intensities.
