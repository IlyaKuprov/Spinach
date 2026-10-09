# examples/kinetics/flux_asymmetric.m

- MATLAB implementation: [examples/kinetics/flux_asymmetric.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/kinetics/flux_asymmetric.m)

- Callable as the no-argument MATLAB function `flux_asymmetric()`; the function constructs its Spinach system and basis before acquisition.

## Purpose and model

This is a two-site intermolecular magnetisation-exchange simulation for two `1H` environments. The spins are separate substances undergoing additive `A+B -> A+B` replacement with swapped atom matching. Concentrations `[2e3 5e2]` and event rate `1` reproduce the former directed rates: 500 from site 1 to site 2 and 2000 in reverse. The initial `L+` state is weighted once by these concentrations; detection is unweighted. Both concentrations are invariant, so the acquisition callback freezes the additive generator at `unit_state` before ordinary linear propagation. No cross-molecular correlations are represented.

## Acquisition and observable

The calculation uses the full `sphten-liouv` basis (`bas.approximation=none`). `liquid` with a frozen-replacement acquisition callback produces an NMR FID, detected with the `1H` `L+` coil. Acquisition uses an empty decoupling list, `offset=900`, `sweep=5000`, 512 points, and zero filling to 1024; the plotted axis is labelled in ppm and inverted. Exponential apodisation uses parameter 6, followed by a shifted FFT. `plot_1d` displays the real spectrum.

## Scope

The source header estimates a calculation time of seconds, not an independently measured runtime. No simulated spectral values or fit result are included in the source; no specific line shape, peak intensity, or fitted flux is claimed.
