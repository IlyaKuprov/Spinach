# examples/kinetics/flux_asymmetric.m

- MATLAB implementation: [examples/kinetics/flux_asymmetric.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/kinetics/flux_asymmetric.m)

- Callable as the no-argument MATLAB function `flux_asymmetric()`; the function constructs its Spinach system and basis before acquisition.

## Purpose and model

This is a two-site, asymmetric intermolecular magnetisation-flux simulation for two `1H` environments. The source sets `sys.magnet=14.1`, scalar Zeeman values `{0,3}`, and directed `inter.chem.flux_rate` entries of 500 from site 1 to site 2 and 2000 from site 2 to site 1, with `inter.chem.flux_type=intermolecular`. It initialises `rho0` with site-1 `L+` weighted by 2000 and site-2 `L+` weighted by 500. The source does not annotate units for the field, scalar values, or flux-rate entries; no units are inferred here.

## Acquisition and observable

The calculation uses the full `sphten-liouv` basis (`bas.approximation=none`). `liquid(spin_system,@acquire,parameters,'nmr')` produces an NMR FID, detected with the `1H` `L+` coil. Acquisition uses an empty decoupling list, `offset=900`, `sweep=5000`, 512 points, and zero filling to 1024; the plotted axis is labelled in ppm and inverted. Exponential apodisation uses parameter 6, followed by a shifted FFT. `plot_1d` displays the real spectrum.

## Scope

The source header estimates a calculation time of seconds, not an independently measured runtime. No simulated spectral values or fit result are included in the source; no specific line shape, peak intensity, or fitted flux is claimed.