# examples/kinetics/exchange_symmetric.m

- MATLAB implementation: [examples/kinetics/exchange_symmetric.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/kinetics/exchange_symmetric.m)

- Callable as the no-argument MATLAB function `exchange_symmetric()`; it creates the Spinach system and basis before running the acquisition.

## Purpose and model

This example simulates a two-site symmetric chemical-exchange NMR pattern for two `1H` environments. The source sets `sys.magnet=14.1`, scalar Zeeman values `{0.0, 3.0}`, exchange sites `{1,2}`, and `inter.chem.rates=[-2e3 2e3; 2e3 -2e3]`, with `inter.chem.concs=[1.0 1.0]`. The field, scalar-value, and rate units are not annotated in this source; the page therefore preserves the configured values without assigning units.

## Acquisition and observable

The full `sphten-liouv` basis is used (`bas.approximation=none`). The initial state is the `1H` chemical-state `L+` operator and the coil is `1H` `L+`. `liquid(spin_system,@acquire,parameters,'nmr')` generates the FID with an empty decoupling list, `offset=900`, `sweep=5000`, 512 points, and zero filling to 1024. The plotted frequency-axis unit is explicitly set to ppm and the axis is inverted. The FID receives exponential apodisation with parameter 6; a shifted FFT is applied and `plot_1d` displays its real part.

## Scope

The source header estimates a calculation time of seconds; this is not a measured runtime. The source defines a simulation and plotting procedure but supplies no numerical spectrum or fitted exchange result, so no line positions, intensities, or fit outcomes are asserted.
