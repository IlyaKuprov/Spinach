# examples/nmr_liquids/sat_rec_strychnine.m

- Signature: `sat_rec_strychnine()`

## Purpose

Simulate a 1H saturation-recovery experiment on strychnine at 250 MHz. Calculation time: minutes.

## Physical / mathematical content

- Loads the strychnine 1H spin system with `strychnine({'1H'})` and sets the magnetic induction to `5.9`.
- Uses Redfield relaxation with `inter.equilibrium='dibari'`, `inter.rlx_keep='kite'`, correlation time `200e-12`, and temperature `298`.
- Uses the `sphten-liouv` formalism with `IK-2` approximation, `scalar_couplings` connectivity, and proximity level `1`.

## Numerical / algorithmic content

- Sets the proximity cutoff to `5.0` and enables `greedy` parallelisation.
- Configures the 1H sequence with offset `1250`, sweep `2500`, `4096` points, maximum delay `0.5`, and `10` delays. The axis is in `ppm` and is inverted.
- Runs `liquid(spin_system,@sat_rec,parameters,'nmr')`, applies exponential apodisation with parameter `6`, then Fourier-transforms and centers the result with `fftshift(fft(fids,[],1))`.
- Plots the real part of the spectra with `plot_1d`.

## Implementation structure

1. Read strychnine spin-system properties and set the magnet, proximity cutoff, and parallelisation option.
2. Set relaxation parameters and the basis; construct the spin system with `create` and `basis`.
3. Set sequence parameters and run the saturation-recovery simulation.
4. Apodise the FIDs, Fourier-transform them, and plot the real spectra.