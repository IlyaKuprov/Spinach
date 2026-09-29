# examples/kinetics/flux_symmetric.m

- MATLAB implementation: [examples/kinetics/flux_symmetric.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/kinetics/flux_symmetric.m)

- Callable as the no-argument MATLAB function `flux_symmetric()`; it constructs the Spinach system and basis before acquisition.

## Purpose and model

This is a two-site symmetric intermolecular magnetisation-flux simulation for two `1H` environments. The source sets `sys.magnet=14.1`, scalar Zeeman values `{0.0,3.0}`, and equal directed `inter.chem.flux_rate` entries of 2000 from site 1 to site 2 and from site 2 to site 1, with `inter.chem.flux_type=intermolecular`. The initial state is the sum of site-1 and site-2 `L+` operators, each weighted by 1.0. Units for the field, scalar values, and flux-rate entries are not stated in this source and are not supplied here.

## Acquisition and observable

The full `sphten-liouv` basis is used (`bas.approximation=none`). The function passes the system, `@acquire` callback, and NMR parameters to `liquid`, with the `1H` `L+` coil and no decoupled spins. Acquisition uses `offset=900`, `sweep=5000`, 512 points, and zero filling to 1024; its plotted axis is labelled in ppm and inverted. The FID is exponentially apodised with parameter 6, transformed with a shifted FFT, and its real spectrum is plotted with `plot_1d`.

## Scope

The source header estimates a calculation time of seconds; this is not a measured runtime. The file specifies a symmetric flux setup rather than a fitted exchange result, and it contains no numerical spectrum from which to report peak positions or intensities.